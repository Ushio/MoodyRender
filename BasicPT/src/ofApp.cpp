#include "alembic_loader.hpp"
#include "render_object.hpp"
#include "microfacet.hpp"
#include "bicubic.hpp"
#include "material.hpp"
#include "value_prportional_sampler.hpp"
#include "geometry.hpp"

#include <xmmintrin.h>
#include <pmmintrin.h>
#include <tbb/tbb.h>
#include <embree3/rtcore.h>

#include "ofApp.h"

#define DEBUG_MODE 0

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


			// Emission
			for (int i = 0; i < _scene->geometries.size(); ++i) {
				const Geometry& g = _scene->geometries[i];
				for (int j = 0; j < g.primitives.size(); ++j) {
					if (g.primitives[j].material->isEmission()) {
						EmissivePrimitiveRef ref;
						ref.geometeryIndex = i;
						ref.primitiveIndex = j;
						auto indices = g.primitives[j].indices;
						ref.area = triArea(g.points[indices.x].P, g.points[indices.y].P, g.points[indices.z].P);
						ref.Ng = triNg(g.points[indices.x].P, g.points[indices.y].P, g.points[indices.z].P, false);
						_emissivePrimitives.push_back(ref);
					}
				}
			}
			_emissionSampler = ValueProportionalSampler<double>(_emissivePrimitives, [](const EmissivePrimitiveRef& ref) { return ref.area; });
		}
		~SceneInterface() {
			rtcReleaseScene(_embreeScene);
			rtcReleaseDevice(_embreeDevice);
		}
		SceneInterface(const SceneInterface &) = delete;
		void operator=(const SceneInterface &) = delete;

		bool occluded(const glm::dvec3 &p, const glm::dvec3 &q, float bias = 1.0e-5f) const {
			glm::dvec3 rd = q - p;
			//double d = glm::distance(p, q);
			//glm::dvec3 rd = (q - p) / d;

			RTCRay ray;
			ray.org_x = p.x;
			ray.org_y = p.y;
			ray.org_z = p.z;
			ray.dir_x = rd.x;
			ray.dir_y = rd.y;
			ray.dir_z = rd.z;
			ray.time = 0.0;

			float tfar = 1.0 - bias;

			ray.tfar = tfar;
			ray.tnear = bias;

			ray.mask = 0;
			ray.id = 0;
			ray.flags = 0;

			rtcOccluded1(_embreeScene, &_context, &ray);
			
			return ray.tfar != tfar;
		}

		bool intersect(const glm::dvec3 &ro, const glm::dvec3 &rd, double tnear, Material *material, double *tmin) const {
			RTCRayHit rayhit;
			rayhit.ray.org_x = ro.x;
			rayhit.ray.org_y = ro.y;
			rayhit.ray.org_z = ro.z;
			rayhit.ray.dir_x = rd.x;
			rayhit.ray.dir_y = rd.y;
			rayhit.ray.dir_z = rd.z;
			rayhit.ray.time = 0.0f;

			rayhit.ray.tfar = (float)*tmin;
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

			(*material)->Ng = glm::normalize(glm::dvec3(rayhit.hit.Ng_x, rayhit.hit.Ng_y, rayhit.hit.Ng_z));

			return true;
		}

		const Camera &camera() const {
			return _scene->camera;
		}

		void sampleEmissiveUniform(PeseudoRandom *random, glm::dvec3 *p, Material *m) const {
			int index = _emissionSampler.sample(random);
			auto prim = _emissivePrimitives[index];
			const Geometry& g = _scene->geometries[prim.geometeryIndex];
			auto material = g.primitives[prim.primitiveIndex].material;
			auto indices = g.primitives[prim.primitiveIndex].indices;

			material->Ng = prim.Ng;
			*m = material;
			*p = uniform_on_triangle(random->uniform(), random->uniform()).evaluate(g.points[indices.x].P, g.points[indices.y].P, g.points[indices.z].P);
		}
		double emissiveArea() const {
			return _emissionSampler.sumValue();
		}

		std::shared_ptr<rt::Scene> _scene;
		RTCDevice _embreeDevice = nullptr;
		RTCScene _embreeScene = nullptr;
		mutable RTCIntersectContext _context;

		std::unordered_map<unsigned int, int> _geomIDToIndex;

		struct EmissivePrimitiveRef {
			int geometeryIndex;
			int primitiveIndex;
			double area;
			glm::dvec3 Ng;
		};
		std::vector<EmissivePrimitiveRef> _emissivePrimitives;
		ValueProportionalSampler<double> _emissionSampler;
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

		void add(int x, int y, glm::dvec3 c) {
			int index = y * _w + x;
			_pixels[index].color += c;
			_pixels[index].sample++;
		}

		struct Pixel {
			int sample = 0;
			glm::dvec3 color;
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

	inline double GTerm(const glm::dvec3 &p, double cosThetaP, const glm::dvec3 &q, double cosThetaQ) {
		return cosThetaP * cosThetaQ / glm::distance2(p, q);
	}

	inline glm::dvec3 radiance(const rt::SceneInterface &scene, glm::dvec3 ro, glm::dvec3 rd, PeseudoRandom *random) {
		glm::dvec3 Lo;
		glm::dvec3 T(1.0);
		for (int i = 0; i < 10; ++i) {
			Material m;
			double tmin = std::numeric_limits<double>::max();
			glm::dvec3 wo = -rd;

			if (scene.intersect(ro, rd, 0.00001, &m, &tmin)) {
				//{
				//	glm::dvec3 p = ro + rd * tmin;

				//	glm::dvec3 q;
				//	Material sampledM;
				//	scene.sampleEmissiveUniform(random, &q, &sampledM);

				//	glm::dvec3 wi = glm::normalize(q - p);
				//	glm::dvec3 emission = bxdf_emission(sampledM, -wi);

				//	glm::dvec3 bxdf = bxdf_evaluate(m, wo, wi);
				//	double pdf_area = (1.0 / scene.emissiveArea());

				//	double cosTheta = glm::dot(bxdf_Ng(m), wi);

				//	double cosThetaP = glm::abs(glm::dot(bxdf_Ng(m), wi));
				//	double cosThetaQ = glm::dot(bxdf_Ng(sampledM), -wi);

				//	double g = GTerm(p, cosThetaP, q, cosThetaQ);

				//	glm::dvec3 contribution = T * bxdf * emission * g;

				//	/* 裏面は発光しない */
				//	if (0.0 < cosThetaQ && glm::any(glm::greaterThanEqual(contribution, glm::dvec3(glm::epsilon<double>())))) {
				//		if (scene.occluded(p, q) == false) {
				//			Lo += contribution / pdf_area;
				//		}
				//	}
				//}
				glm::dvec3 wi = m->sample(random, wo);
				glm::dvec3 bxdf = m->bxdf(wo, wi);
				glm::dvec3 emission = m->emission(wo);
				double pdf = m->pdf(wo, wi);
				double cosTheta = std::abs(glm::dot(m->Ng, wi));

				// MicrofacetConductorMaterial の ダメサンプルを数える
				if (m.get<MicrofacetConductorMaterial>()) {
					static std::atomic<int> badSample = 0;
					static std::atomic<int> numSample = 0;
					if (glm::dot(m->Ng, wi) < 0.0) {
						badSample++;
					}
					numSample++;
					if (numSample % 100000 == 0) {
						printf("bad sample %.2f%%\n", (double)badSample / numSample * 100.0);
					}
				}

				//if (i == 0) {
				//	Lo += emission * T;
				//}
				Lo += emission * T;

				if (glm::any(glm::greaterThanEqual(bxdf, glm::dvec3(1.0e-6f)))) {
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

	//inline glm::dvec3 radiance(const rt::SceneInterface &scene, glm::dvec3 ro, glm::dvec3 rd, PeseudoRandom *random) {
	//	glm::dvec3 Lo;
	//	glm::dvec3 T(1.0);
	//	for (int i = 0; i < 10; ++i) {
	//		Material m;
	//		double tmin = std::numeric_limits<double>::max();
	//		glm::dvec3 wo = -rd;

	//		if (scene.intersect(ro, rd, 0.00001, &m, &tmin)) {
	//			glm::dvec3 wi = bxdf_sample(m, random, wo);
	//			glm::dvec3 bxdf = bxdf_evaluate(m, wo, wi);
	//			glm::dvec3 emission = bxdf_emission(m, wo);
	//			double pdf = bxdf_pdf(m, wo, wi);
	//			double cosTheta = glm::dot(bxdf_Ng(m), wi);

	//			Lo += emission * T;

	//			if(glm::any(glm::greaterThanEqual(bxdf, glm::dvec3(1.0e-6f)))) {
	//				T *= bxdf * cosTheta / pdf;
	//			}
	//			else {
	//				break;
	//			}

	//			ro = (ro + rd * tmin);
	//			rd = wi;
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

#if DEBUG_MODE
			int focusX = 200;
			int focusY = 200;

			for (int y = 0; y < _scene->camera.imageHeight(); ++y) {
				for (int x = 0; x < _scene->camera.imageWidth(); ++x) {
					if (x != focusX || y != focusY) {
						continue;
					}
					// PeseudoRandom *random = &_image.pixel(x, y)->random;
					auto cp = _image.pixel(x, y)->random;
					PeseudoRandom *random = &cp;

					glm::dvec3 o;
					glm::dvec3 d;
					_scene->camera.sampleRay(random, x, y, &o, &d);

					auto r = radiance(*_sceneInterface, o, d, random);

					for (int i = 0; i < r.length(); ++i) {
						if (glm::isfinite(r[i]) == false) {
							r[i] = 0.0;
						}
						if (r[i] < 0.0 || 1000.0 < r[i]) {
							r[i] = 0.0;
						}
					}
					_image.add(x, y, r);
				}
			}
#else
			tbb::parallel_for(tbb::blocked_range<int>(0, _scene->camera.imageHeight()), [&](const tbb::blocked_range<int> &range) {
				for (int y = range.begin(); y < range.end(); ++y) {
					for (int x = 0; x < _scene->camera.imageWidth(); ++x) {
						PeseudoRandom *random = &_image.pixel(x, y)->random;
						glm::dvec3 o;
						glm::dvec3 d;
						_scene->camera.sampleRay(random, x, y, &o, &d);

						auto r = radiance(*_sceneInterface, o, d, random);

						for (int i = 0; i < r.length(); ++i) {
							if (glm::isfinite(r[i]) == false) {
								r[i] = 0.0;
							}
							if (r[i] < 0.0 || 1000.0 < r[i]) {
								r[i] = 0.0;
							}
						}
						_image.add(x, y, r);
					}
				}
			});
#endif
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

	double scale = 1.0;
	for (int y = 0; y < image.height(); ++y) {
		for (int x = 0; x < image.width(); ++x) {
			int index = y * image.width() + x;
			const auto &px = *image.pixel(x, y);
			auto L = px.color / (double)px.sample;
			dst[index * 3 + 0] = (uint8_t)glm::clamp(glm::pow(L.x * scale, 1.0 / 2.2) * 255.0, 0.0, 255.99999);
			dst[index * 3 + 1] = (uint8_t)glm::clamp(glm::pow(L.y * scale, 1.0 / 2.2) * 255.0, 0.0, 255.99999);
			dst[index * 3 + 2] = (uint8_t)glm::clamp(glm::pow(L.z * scale, 1.0 / 2.2) * 255.0, 0.0, 255.99999);
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

	_camera.setNearClip(0.1);
	_camera.setFarClip(100.0);
	_camera.setDistance(5.0);

	scene = std::shared_ptr<rt::Scene>(new rt::Scene());
	rt::loadFromABC(ofToDataPath("cornelbox.abc").c_str(), *scene);
	// rt::loadFromABC(ofToDataPath("cornelbox.abc").c_str(), *scene);

	renderer = std::shared_ptr<rt::PTRenderer>(new rt::PTRenderer(scene));

	rt::CoupledBRDFConductor::load(ofToDataPath("baked/albedo_specular_conductor.xml").c_str(), ofToDataPath("baked/albedo_specular_conductor_avg.xml").c_str());
	rt::CoupledBRDFDielectrics::load(ofToDataPath("baked/albedo_specular_dielectrics.xml").c_str(), ofToDataPath("baked/albedo_specular_dielectrics_avg.xml").c_str());
}

//--------------------------------------------------------------
void ofApp::update() {

}

//--------------------------------------------------------------
void ofApp::draw() {
	ofEnableDepthTest();

	ofClear(0);
	_camera.begin();
	ofPushMatrix();
	ofRotateZ(90.0);
	ofSetColor(255);
	ofDrawGridPlane(1.0);
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
		ofDrawSphere(origin.x, origin.y, origin.z, 0.05);

		ofSetColor(0, 0, 255);
		auto frontP = origin + camera.front() * 0.2;
		ofDrawLine(origin.x, origin.y, origin.z, frontP.x, frontP.y, frontP.z);

		ofSetColor(255, 0, 0);
		auto upP = origin + camera.up() * 0.2;
		ofDrawLine(origin.x, origin.y, origin.z, upP.x, upP.y, upP.z);
	}

	{
		rt::Xor64 random;
		int x = ofMap(ofGetMouseX(), 0, ofGetWidth(), 0, scene->camera.imageWidth());
		int y = ofMap(ofGetMouseY(), 0, ofGetHeight(), 0, scene->camera.imageHeight());
		glm::dvec3 o;
		glm::dvec3 d;
		scene->camera.sampleRay(&random, x, y, &o, &d);

		rt::Material m;
		double tmin = std::numeric_limits<double>::max();
		if (renderer->sceneInterface().intersect(o, d, 0.0, &m, &tmin)) {
			ofSetColor(255, 0, 0);
			auto p = o + d * tmin;
			ofDrawLine(o.x, o.y, o.z, p.x, p.y, p.z);
			
			auto pn = p + m->Ng * 0.1;
			ofDrawLine(p.x, p.y, p.z, pn.x, pn.y, pn.z);
		} else {
			ofSetColor(255);
			auto p = o + d * 10.0;
			ofDrawLine(o.x, o.y, o.z, p.x, p.y, p.z);
		}
	}

	{
		renderer->step();

		if (ofGetFrameNum() % 5 == 0) {
			_image.setFromPixels(toOf(renderer->_image));
		}
		if (renderer->stepCount() == 128) {
			_image.setFromPixels(toOf(renderer->_image));
			_image.save("128spp.png");
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

	//{
	//	static rt::Xor64 random;


	//	ofMesh mesh;
	//	mesh.setMode(OF_PRIMITIVE_POINTS);

	//	ofSetColor(255, 0, 0);

	//	for (int i = 0; i < 3000; ++i) {
	//		glm::dvec3 p;
	//		rt::Material m;
	//		renderer->sceneInterface().sampleEmissiveUniform(&random, &p, &m);
	//		mesh.addVertex(p);

	//		ofDrawLine(p, p + rt::bxdf_Ng(m) * 0.2);
	//	}

	//	mesh.draw();
	//}

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
		_camera.setNearClip(0.1);
		_camera.setFarClip(100.0);
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
