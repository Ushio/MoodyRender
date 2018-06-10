﻿#include "alembic_loader.hpp"
#include "render_object.hpp"
#include "microfacet.hpp"
#include "bicubic.hpp"
#include "material.hpp"

#include <xmmintrin.h>
#include <pmmintrin.h>
#include <tbb/tbb.h>
#include <embree3/rtcore.h>

#include "ofApp.h"

namespace rt {
	inline void EmbreeErorrHandler(void* userPtr, RTCError code, const char* str) {
		printf("Embree Error [%d] %s\n", code, str);
	}
	class SceneInterface {
	public:
		SceneInterface(std::shared_ptr<rt::Scene> scene):_scene(scene) {
			_embreeDevice = rtcNewDevice(nullptr);
			rtcSetDeviceErrorFunction(_embreeDevice, EmbreeErorrHandler, nullptr);
			_embreeScene = rtcNewScene(_embreeDevice);
			rtcSetSceneBuildQuality(_embreeScene, RTC_BUILD_QUALITY_HIGH);
			// RTC_BUILD_QUALITY_LOW, RTC_BUILD_QUALITY_MEDIUM, RTC_BUILD_QUALITY_HIGH

			for (int i = 0; i < _scene->geometries.size(); ++i) {
				RTCGeometry embreeGeometry = rtcNewGeometry(_embreeDevice, RTC_GEOMETRY_TYPE_TRIANGLE);

				// バッファーはGeometryに結びつき、所有される
				float *pVertexBuffer = (float *)rtcSetNewGeometryBuffer(embreeGeometry, RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT3, sizeof(float) * 3, _scene->geometries[i].points.size());
				for (int j = 0; j < _scene->geometries[i].points.size(); ++j) {
					auto p = _scene->geometries[i].points[j].P;

					for (int k = 0; k < 3; ++k) {
						pVertexBuffer[j * 3 + k] = p[k];
					}
				}
				uint32_t *pIndexBuffer = (uint32_t *)rtcSetNewGeometryBuffer(embreeGeometry, RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT3, sizeof(uint32_t) * 3, _scene->geometries[i].primitives.size());
				for (int j = 0; j < _scene->geometries[i].primitives.size(); ++j) {
					auto indices = _scene->geometries[i].primitives[j].indices;
					for (int k = 0; k < 3; ++k) {
						pIndexBuffer[j * 3 + k] = indices[k];
					}
				}
				rtcCommitGeometry(embreeGeometry);
				auto geomID = rtcAttachGeometry(_embreeScene, embreeGeometry);
				_geomIDToIndex[geomID] = i;
				rtcReleaseGeometry(embreeGeometry);
			}

			rtcCommitScene(_embreeScene);

			rtcInitIntersectContext(&_context);
		}
		~SceneInterface() {
			rtcReleaseScene(_embreeScene);
			rtcReleaseDevice(_embreeDevice);
		}
		SceneInterface(const SceneInterface &) = delete;
		void operator=(const SceneInterface &) = delete;

		bool intersect(const glm::vec3 &ro, const glm::vec3 &rd, float tnear, Material *material, float *tmin) const {
			RTCRayHit rayhit;
			rayhit.ray.org_x = ro.x;
			rayhit.ray.org_y = ro.y;
			rayhit.ray.org_z = ro.z;
			rayhit.ray.dir_x = rd.x;
			rayhit.ray.dir_y = rd.y;
			rayhit.ray.dir_z = rd.z;
			rayhit.ray.time = 0.0f;

			rayhit.ray.tfar = *tmin;
			rayhit.ray.tnear = tnear;

			rayhit.ray.mask = 0;
			rayhit.ray.id = 0;
			rayhit.ray.flags = 0;
			rayhit.hit.geomID = RTC_INVALID_GEOMETRY_ID;
			rayhit.hit.instID[0] = RTC_INVALID_GEOMETRY_ID;
			rtcIntersect1(_embreeScene, &_context, &rayhit);

			if (rayhit.hit.geomID == RTC_INVALID_GEOMETRY_ID) {
				return false;
			}

			*tmin = rayhit.ray.tfar;
			
			int index = _geomIDToIndex.find(rayhit.hit.geomID)->second;
			const auto &prim = _scene->geometries[index].primitives[rayhit.hit.primID];
			*material = prim.material;
			glm::vec3 NgUnnormalized(rayhit.hit.Ng_x, rayhit.hit.Ng_y, rayhit.hit.Ng_z);

			strict_variant::apply_visitor(MaterialVisitor::SetNg(glm::normalize(NgUnnormalized)), *material);

			return true;
		}

		const Camera &camera() const {
			return _scene->camera;
		}

		std::shared_ptr<rt::Scene> _scene;
		RTCDevice _embreeDevice = nullptr;
		RTCScene _embreeScene = nullptr;
		mutable RTCIntersectContext _context;

		std::unordered_map<unsigned int, int> _geomIDToIndex;
	};

	class Image {
	public:
		Image(int w, int h) :_w(w), _h(h), _pixels(h*w) {
			for (int i = 0; i < _pixels.size(); ++i) {
				_pixels[i].random = Xor64(i + 1);
				for (int j = 0; j < 100; ++j) {
					_pixels[i].random.generate();
				}
			}
		}
		int width() const {
			return _w;
		}
		int height() const {
			return _h;
		}

		void add(int x, int y, glm::vec3 c) {
			int index = y * _w + x;
			_pixels[index].color += c;
			_pixels[index].sample++;
		}

		struct Pixel {
			int sample = 0;
			glm::vec3 color;
			Xor64 random;
		};
		const Pixel *pixel(int x, int y) const {
			return _pixels.data() + y * _w + x;
		}
		Pixel *pixel(int x, int y) {
			return _pixels.data() + y * _w + x;
		}
	private:
		int _w = 0;
		int _h = 0;
		std::vector<Pixel> _pixels;
	};

	inline glm::vec3 radiance(const rt::SceneInterface &scene, glm::vec3 ro, glm::vec3 rd, PeseudoRandom *random) {
		glm::vec3 Lo;
		glm::vec3 T(1.0);
		for (int i = 0; i < 10; ++i) {
			Material m;
			float tmin = std::numeric_limits<float>::max();
			glm::vec3 wo = -rd;

			if (scene.intersect(ro, rd, 0.00001f, &m, &tmin)) {
				glm::vec3 wi = bxdf_sample(m, random, wo);
				glm::vec3 bxdf = bxdf_evaluate(m, wo, wi);
				glm::vec3 emission = bxdf_emission(m, wo);
				float pdf = bxdf_pdf(m, wo, wi);
				float cosTheta = glm::dot(bxdf_Ng(m), wi);

				Lo += emission * T;

				if(glm::any(glm::greaterThanEqual(bxdf, glm::vec3(1.0e-6f)))) {
					T *= bxdf * cosTheta / pdf;
				}
				else {
					break;
				}

				ro = (ro + rd * tmin);
				rd = wi;
			}
			else {
				break;
			}
		}
		return Lo;
	}

	//inline glm::vec3 radiance(const rt::SceneInterface &scene, glm::vec3 ro, glm::vec3 rd, PeseudoRandom *random) {
	//	glm::vec3 Lo;
	//	glm::vec3 T(1.0);
	//	for (int i = 0; i < 10; ++i) {
	//		Material mat;
	//		float tmin = std::numeric_limits<float>::max();

	//		if (scene.intersect(ro, rd, 0.00001f, &mat, &tmin)) {
	//			if (auto material = strict_variant::get<LambertianMaterial>(&mat)) {
	//				glm::vec3 wi = LambertianSampler::sample(random, material->Ng);
	//				float cos_term = glm::dot(wi, material->Ng);

	//				//glm::vec3 sample = sample_cosine_weighted_hemisphere_brdf(random);
	//				//float pdf_omega = cosine_weighted_hemisphere_pdf_brdf(sample);
	//				//glm::vec3 wi = from_bxdf(material->Ng, sample);

	//				glm::vec3 brdf = material->R * glm::vec3(glm::one_over_pi<float>());
	//				//float cos_term = abs_cos_theta_bxdf(sample);

	//				Lo += material->Le * T;

	//				T *= brdf * cos_term / LambertianSampler::pdf(wi, material->Ng);

	//				ro = (ro + rd * tmin);
	//				rd = wi;
	//			}
	//			else if (auto material = strict_variant::get<MicrofacetCoupledConductorMaterial>(&mat)) {
	//				glm::vec3 wo = -rd;
	//				float alpha = 0.3f;

	//				/*
	//				 コサイン重点サンプリング
	//				*/
	//				//glm::vec3 sample = sample_cosine_weighted_hemisphere_brdf(random);
	//				//float pdf_omega = cosine_weighted_hemisphere_pdf_brdf(sample);
	//				//glm::vec3 wi = from_bxdf(material->Ng, sample);

	//				/*
	//				 ハーフベクトルの重点サンプリング
	//				*/
	//				//float theta = std::atan(std::sqrt(-alpha * alpha * std::log(1.0f - random->uniform())));
	//				//float phi = random->uniform(0.0f, glm::two_pi<double>());
	//				//glm::vec3 sample = polar_to_cartesian(theta, phi);
	//				//glm::vec3 harf = from_bxdf(material->Ng, sample);
	//				//glm::vec3 wi = glm::reflect(-wo, harf);
	//				//float pdf_omega = D_Beckmann(material->Ng, harf, alpha) * glm::dot(material->Ng, harf) / (4.0f * glm::dot(wi, harf));
	//				//if (glm::dot(material->Ng, wi) <= 0.0f) {
	//				//	T = glm::vec3(0.0);
	//				//	break;
	//				//}

	//				// 
	//				//glm::vec3 wi = NDFImportanceSampler::sample_wi_Beckmann(random, alpha, wo, material->Ng);
	//				//float pdf_omega = NDFImportanceSampler::pdfBeckmann(wi, alpha, wo, material->Ng);

	//				// ミックス 重点サンプリング
	//				glm::vec3 wi;
	//				// float spAlbedo = CoupledBRDFDielectrics::specularAlbedo().sample(alpha, glm::dot(material->Ng, wo));
	//				float spAlbedo = CoupledBRDFConductor::specularAlbedo().sample(alpha, glm::dot(material->Ng, wo));

	//				if (random->uniformf() < spAlbedo) {
	//					wi = BeckmannImportanceSampler::sample(random, alpha, wo, material->Ng);
	//				}
	//				else {
	//					wi = LambertianSampler::sample(random, material->Ng);
	//				}

	//				float pdf_omega = 
	//					spAlbedo * BeckmannImportanceSampler::pdf(wi, alpha, wo, material->Ng)
	//					+
	//					(1.0f - spAlbedo) * LambertianSampler::pdf(wi, material->Ng);

	//				glm::vec3 h = glm::normalize(wi + wo);
	//				float d = D_Beckmann(material->Ng, h, alpha);
	//				float g = G2_height_correlated_beckmann(wi, wo, h, material->Ng, alpha);

	//				float cos_term_wo = glm::dot(material->Ng, wo);
	//				float cos_term_wi = glm::dot(material->Ng, wi);

	//				float brdf_without_f = d * g / (4.0f * cos_term_wo * cos_term_wi);

	//				glm::vec3 eta(0.15557f, 0.42415f, 1.3821f);
	//				glm::vec3 k(3.6024f, 2.4721f, 1.9155f);

	//				float cosThetaFresnel = glm::dot(h, wo);
	//				glm::vec3 f = glm::vec3(
	//					fresnel_unpolarized(eta.r, k.r, cosThetaFresnel),
	//					fresnel_unpolarized(eta.g, k.g, cosThetaFresnel),
	//					fresnel_unpolarized(eta.b, k.b, cosThetaFresnel)
	//				);
	//				// glm::vec3 f = glm::vec3(fresnel_dielectrics(cosThetaFresnel));

	//				glm::vec3 brdf_spec = f * brdf_without_f;

	//				// glm::vec3 kLambda(1.0f, 0.447067, 0.246); // 適当

	//				glm::vec3 kLambda(
	//					fresnel_unpolarized(eta.r, k.r, 0.0),
	//					fresnel_unpolarized(eta.g, k.g, 0.0),
	//					fresnel_unpolarized(eta.b, k.b, 0.0)
	//				);

	//				//glm::vec3 brdf_diff = kLambda
	//				//	* (1.0f - CoupledBRDFDielectrics::specularAlbedo().sample(alpha, cos_term_wo))
	//				//	* (1.0f - CoupledBRDFDielectrics::specularAlbedo().sample(alpha, cos_term_wi))
	//				//	/ (glm::pi<float>() * (1.0f - CoupledBRDFDielectrics::specularAvgAlbedo().sample(alpha)));
	//				glm::vec3 brdf_diff = kLambda
	//					* (1.0f - CoupledBRDFConductor::specularAlbedo().sample(alpha, cos_term_wo))
	//					* (1.0f - CoupledBRDFConductor::specularAlbedo().sample(alpha, cos_term_wi))
	//					/ (glm::pi<float>() * (1.0f - CoupledBRDFConductor::specularAvgAlbedo().sample(alpha)));

	//				// glm::vec3 brdf = brdf_spec + brdf_diff;
	//				glm::vec3 brdf = brdf_spec;

	//				if (glm::dot(material->Ng, wi) <= 0.0f) {
	//					brdf = glm::vec3();
	//				}

	//				T *= brdf * cos_term_wi / pdf_omega;

	//				ro = (ro + rd * tmin);
	//				rd = wi;
	//			}
	//			else if (auto material = strict_variant::get<SpecularMaterial>(&mat)) {
	//				glm::vec3 wi = glm::reflect(rd, material->Ng);
	//				ro = (ro + rd * tmin);
	//				rd = wi;
	//			}
	//		}
	//		else {
	//			break;
	//		}
	//	}
	//	return Lo;
	//}

	class PTRenderer {
	public:
		PTRenderer(std::shared_ptr<rt::Scene> scene)
			: _scene(scene)
			, _sceneInterface(new rt::SceneInterface(scene))
			, _image(scene->camera.imageWidth(), scene->camera.imageHeight()) {

		}
		void step() {
			_steps++;

			tbb::parallel_for(tbb::blocked_range<int>(0, _scene->camera.imageHeight()), [&](const tbb::blocked_range<int> &range) {
				for (int y = range.begin(); y < range.end(); ++y) {
					for (int x = 0; x < _scene->camera.imageWidth(); ++x) {
						PeseudoRandom *random = &_image.pixel(x, y)->random;
						glm::vec3 o;
						glm::vec3 d;
						_scene->camera.sampleRay(random, x, y, &o, &d);

						auto r = radiance(*_sceneInterface, o, d, random);

						for (int i = 0; i < r.length(); ++i) {
							if (glm::isfinite(r[i]) == false) {
								r[i] = 0.0f;
							}
							if (r[i] < 0.0f || 1000.0f < r[i]) {
								r[i] = 0.0f;
							}
						}
						_image.add(x, y, r);
					}
				}
			});
		}
		int stepCount() const {
			return _steps;
		}

		const rt::SceneInterface &sceneInterface() const {
			return *_sceneInterface;
		}

		std::shared_ptr<rt::Scene> _scene;
		std::shared_ptr<rt::SceneInterface> _sceneInterface;
		Image _image;
		int _steps = 0;
	};
}

// PT
inline ofPixels toOf(const rt::Image &image) {
	ofPixels pixels;
	pixels.allocate(image.width(), image.height(), OF_IMAGE_COLOR);
	uint8_t *dst = pixels.getPixels();

	float scale = 1.0f;
	for (int y = 0; y < image.height(); ++y) {
		for (int x = 0; x < image.width(); ++x) {
			int index = y * image.width() + x;
			const auto &px = *image.pixel(x, y);
			auto L = px.color / (float)px.sample;
			dst[index * 3 + 0] = (uint8_t)glm::clamp(glm::pow(L.x * scale, 1.0f / 2.2f) * 255.0f, 0.0f, 255.99999f);
			dst[index * 3 + 1] = (uint8_t)glm::clamp(glm::pow(L.y * scale, 1.0f / 2.2f) * 255.0f, 0.0f, 255.99999f);
			dst[index * 3 + 2] = (uint8_t)glm::clamp(glm::pow(L.z * scale, 1.0f / 2.2f) * 255.0f, 0.0f, 255.99999f);
		}
	}
	return pixels;
}


std::shared_ptr<rt::Scene> scene;
std::shared_ptr<rt::PTRenderer> renderer;

//--------------------------------------------------------------
void ofApp::setup(){
	

	_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
	_MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);

	ofSetVerticalSync(false);

	_camera.setNearClip(0.1f);
	_camera.setFarClip(100.0f);
	_camera.setDistance(5.0f);

	scene = std::shared_ptr<rt::Scene>(new rt::Scene());
	rt::loadFromABC(ofToDataPath("cornelbox.abc").c_str(), *scene);
	// rt::loadFromABC(ofToDataPath("cornelbox.abc").c_str(), *scene);

	renderer = std::shared_ptr<rt::PTRenderer>(new rt::PTRenderer(scene));

	rt::CoupledBRDFConductor::load(ofToDataPath("baked/albedo_specular_conductor.xml").c_str(), ofToDataPath("baked/albedo_specular_conductor_avg.xml").c_str());
	rt::CoupledBRDFDielectrics::load(ofToDataPath("baked/albedo_specular_dielectrics.xml").c_str(), ofToDataPath("baked/albedo_specular_dielectrics_avg.xml").c_str());
}

//--------------------------------------------------------------
void ofApp::update() {
	renderer->step();

	if (ofGetFrameNum() % 5 == 0) {
		_image.setFromPixels(toOf(renderer->_image));
	}

	if (renderer->stepCount() == 512) {
		_image.setFromPixels(toOf(renderer->_image));
		_image.save("512spp.png");
	}
	if (renderer->stepCount() == 1024) {
		_image.setFromPixels(toOf(renderer->_image));
		_image.save("1024spp.png");
	}
}

//--------------------------------------------------------------
void ofApp::draw() {
	ofEnableDepthTest();

	ofClear(0);
	_camera.begin();
	ofPushMatrix();
	ofRotateZ(90.0f);
	ofSetColor(255);
	ofDrawGridPlane(1.0f);
	ofPopMatrix();

	ofPushMatrix();
	ofDrawAxis(50);
	ofPopMatrix();

	ofSetColor(255);
	for (int i = 0; i < scene->geometries.size(); ++i) {
		ofMesh mesh;

		auto geometry = scene->geometries[i];
		for (int j = 0; j < geometry.points.size(); ++j) {
			auto p = geometry.points[j].P;
			mesh.addVertex(ofVec3f(p.x, p.y, p.z));
		}
		for (int j = 0; j < geometry.primitives.size(); ++j) {
			auto prim = geometry.primitives[j];
			mesh.addIndex(prim.indices[0]);
			mesh.addIndex(prim.indices[1]);
			mesh.addIndex(prim.indices[2]);
		}
		mesh.drawWireframe();
	}
	{
		auto camera = scene->camera;
		ofSetColor(255);
		auto origin = camera.origin();
		ofDrawSphere(origin.x, origin.y, origin.z, 0.05f);

		ofSetColor(0, 0, 255);
		auto frontP = origin + camera.front() * 0.2f;
		ofDrawLine(origin.x, origin.y, origin.z, frontP.x, frontP.y, frontP.z);

		ofSetColor(255, 0, 0);
		auto upP = origin + camera.up() * 0.2f;
		ofDrawLine(origin.x, origin.y, origin.z, upP.x, upP.y, upP.z);
	}

	{
		rt::Xor64 random;
		int x = ofMap(ofGetMouseX(), 0, ofGetWidth(), 0, scene->camera.imageWidth());
		int y = ofMap(ofGetMouseY(), 0, ofGetHeight(), 0, scene->camera.imageHeight());
		glm::vec3 o;
		glm::vec3 d;
		scene->camera.sampleRay(&random, x, y, &o, &d);

		rt::Material m;
		float tmin = std::numeric_limits<float>::max();
		if (renderer->sceneInterface().intersect(o, d, 0.0f, &m, &tmin)) {
			ofSetColor(255, 0, 0);
			auto p = o + d * tmin;
			ofDrawLine(o.x, o.y, o.z, p.x, p.y, p.z);
			
			auto pn = p + strict_variant::apply_visitor(rt::MaterialVisitor::GetNg(), m) * 0.1f;
			ofDrawLine(p.x, p.y, p.z, pn.x, pn.y, pn.z);
		} else {
			ofSetColor(255);
			auto p = o + d * 10.0f;
			ofDrawLine(o.x, o.y, o.z, p.x, p.y, p.z);
		}
	}

	_camera.end();

	ofDisableDepthTest();
	ofSetColor(255);


	if (_image.isAllocated() && _showImage) {
		_image.draw(0, 0);
	}

	char buffer[256];
	sprintf(buffer, "%d sample, fps = %.3f", renderer->stepCount(), ofGetFrameRate());
	ofDrawBitmapString(buffer, 10, 10);
}

//--------------------------------------------------------------
void ofApp::keyPressed(int key) {
	if (key == 'f') {
		auto camera = scene->camera;
		ofSetWindowShape(camera.imageWidth(), camera.imageHeight());
		_camera.setNearClip(0.1f);
		_camera.setFarClip(100.0f);
		_camera.setFov(glm::degrees(camera.setting().fovy));
		_camera.setPosition(camera.origin().x, camera.origin().y, camera.origin().z);

		auto lookAt = camera.origin() + camera.front();
		_camera.lookAt(ofVec3f(lookAt.x, lookAt.y, lookAt.z), ofVec3f(camera.up().x, camera.up().y, camera.up().z));
	}

	if (key == ' ') {
		_showImage = !_showImage;
	}

	if (key == 's') {
		_image.save("pt.png");
		ofFloatImage img(_image);
	}
}

//--------------------------------------------------------------
void ofApp::keyReleased(int key){

}

//--------------------------------------------------------------
void ofApp::mouseMoved(int x, int y ){

}

//--------------------------------------------------------------
void ofApp::mouseDragged(int x, int y, int button){

}

//--------------------------------------------------------------
void ofApp::mousePressed(int x, int y, int button){

}

//--------------------------------------------------------------
void ofApp::mouseReleased(int x, int y, int button){

}

//--------------------------------------------------------------
void ofApp::mouseEntered(int x, int y){

}

//--------------------------------------------------------------
void ofApp::mouseExited(int x, int y){

}

//--------------------------------------------------------------
void ofApp::windowResized(int w, int h){

}

//--------------------------------------------------------------
void ofApp::gotMessage(ofMessage msg){

}

//--------------------------------------------------------------
void ofApp::dragEvent(ofDragInfo dragInfo){ 

}
