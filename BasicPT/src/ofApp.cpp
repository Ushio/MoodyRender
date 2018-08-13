#include "alembic_loader.hpp"
#include "render_object.hpp"
#include "microfacet.hpp"
#include "bicubic.hpp"
#include "material.hpp"
#include "value_prportional_sampler.hpp"
#include "geometry.hpp"
#include "stopwatch.hpp"
#include "direct_sampler.hpp"

#include <random>
#include <xmmintrin.h>
#include <pmmintrin.h>
#include <tbb/tbb.h>
#include <embree3/rtcore.h>

#include "ofApp.h"
#include "ofxImGuiLite.hpp"

#define DEBUG_MODE 0

namespace rt {
	
}

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
				rtcAttachGeometryByID(_embreeScene, embreeGeometry, i);
				rtcReleaseGeometry(embreeGeometry);
			}

			rtcCommitScene(_embreeScene);

			rtcInitIntersectContext(&_context);


			for (int i = 0; i < _scene->geometries.size(); ++i) {
				Geometry& g = _scene->geometries[i];
				SphericalRectangleSampler *previous_sr_sampler = nullptr;

				for (int j = 0; j < g.primitives.size(); ++j) {
					Geometry::Primitive &p = g.primitives[j];
					IMaterial *m = g.primitives[j].material.get();
					LambertianMaterial *lambertian = dynamic_cast<LambertianMaterial *>(m);
					if (lambertian && lambertian->isEmission()) {
						if (auto sample = lambertian->samplingStrategy.get<AreaSample>()) {
							auto a = g.points[p.indices[0]].P;
							auto b = g.points[p.indices[1]].P;
							auto c = g.points[p.indices[2]].P;
							TriangleAreaSampler *sampler = new TriangleAreaSampler(a, b, c, lambertian->backEmission, lambertian->Le);
							lambertian->sampler = sampler;
							_directSamplers.emplace_back(sampler);
						}
						else if (auto sample = lambertian->samplingStrategy.get<SphericalRectangleSample>())
						{
							// あまり綺麗ではないが、triangleIndex 0, 1, 0, 1...となる決まりにする。
							if (sample->triangleIndex == 0) {
								SphericalRectangleSampler *sampler = new SphericalRectangleSampler(sample->s, sample->ex, sample->ey, lambertian->backEmission, lambertian->Le);
								previous_sr_sampler = sampler;
								lambertian->sampler = sampler;
								_directSamplers.emplace_back(sampler);
							}
							else {
								lambertian->sampler = previous_sr_sampler;
							}
						}
					}
				}
			}

			RTCBounds bounds;
			rtcGetSceneBounds(_embreeScene, &bounds);

			float maxWide = bounds.upper_x - bounds.lower_x;
			maxWide = std::max(maxWide, bounds.upper_y - bounds.lower_y);
			maxWide = std::max(maxWide, bounds.upper_z - bounds.lower_z);

			_sceneAdaptiveEps = std::max(maxWide * 1.0e-5, 1.0e-5);
			// printf("_sceneAdaptiveEps %.10f\n", _sceneAdaptiveEps);
		}
		~SceneInterface() {
			rtcReleaseScene(_embreeScene);
			rtcReleaseDevice(_embreeDevice);
		}
		SceneInterface(const SceneInterface &) = delete;
		void operator=(const SceneInterface &) = delete;

		double adaptiveEps() const {
			return _sceneAdaptiveEps;
		}

		bool occluded(const glm::dvec3 &p, const glm::dvec3 &q) const {
			glm::dvec3 rd = q - p;

			RTCRay ray;
			ray.org_x = p.x;
			ray.org_y = p.y;
			ray.org_z = p.z;
			ray.dir_x = rd.x;
			ray.dir_y = rd.y;
			ray.dir_z = rd.z;
			ray.time = 0.0;

			ray.tfar = 1.0f;
			ray.tnear = 0.0f;

			ray.mask = 0;
			ray.id = 0;
			ray.flags = 0;

			rtcOccluded1(_embreeScene, &_context, &ray);
			
			return ray.tfar != 1.0f;
		}

		bool intersect(const glm::dvec3 &ro, const glm::dvec3 &rd, Material *material, float *tmin) const {
			RTCRayHit rayhit;
			rayhit.ray.org_x = ro.x;
			rayhit.ray.org_y = ro.y;
			rayhit.ray.org_z = ro.z;
			rayhit.ray.dir_x = rd.x;
			rayhit.ray.dir_y = rd.y;
			rayhit.ray.dir_z = rd.z;
			rayhit.ray.time = 0.0f;

			rayhit.ray.tfar = FLT_MAX;
			rayhit.ray.tnear = 0.0f;

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
			
			int index = rayhit.hit.geomID;
			const auto &geom = _scene->geometries[index];
			const auto &prim = geom.primitives[rayhit.hit.primID];
			*material = prim.material;

			// TODO これは後で削る
			glm::dvec3 Ng = glm::normalize(glm::dvec3(rayhit.hit.Ng_x, rayhit.hit.Ng_y, rayhit.hit.Ng_z));

			// 裏面
			bool backfacing = false;
			if (glm::dot(rd, Ng) > 0.0)
			{
				Ng = -Ng;
				backfacing = true;
			}

			(*material)->Ng = Ng;
			(*material)->backfacing = backfacing;

			/*
			https://embree.github.io/api.html
			t_uv = (1-u-v)*t0 + u*t1 + v*t2
				 = t0 + u*(t1-t0) + v*(t2-t0)
			*/
			//double u = rayhit.hit.u;
			//double v = rayhit.hit.v;
			//auto v0 = geom.points[prim.indices[0]].P;
			//auto v1 = geom.points[prim.indices[1]].P;
			//auto v2 = geom.points[prim.indices[2]].P;
			//(*material)->p = (1.0 - u - v) * v0 + u * v1 + v * v2;

			(*material)->p = ro + rd * double(*tmin);
			return true;
		}

		const Camera &camera() const {
			return _scene->camera;
		}

		std::vector<std::unique_ptr<IDirectSampler>>::const_iterator sampler_begin() const {
			return _directSamplers.begin();
		}
		std::vector<std::unique_ptr<IDirectSampler>>::const_iterator sampler_end() const {
			return _directSamplers.end();
		}
		int samplerCount() const {
			return _directSamplers.size();
		}

		std::shared_ptr<rt::Scene> _scene;
		RTCDevice _embreeDevice = nullptr;
		RTCScene _embreeScene = nullptr;
		mutable RTCIntersectContext _context;

		std::vector<std::unique_ptr<IDirectSampler>> _directSamplers;

		double _sceneAdaptiveEps = 0.0f;
	};

	class Image {
	public:
		Image(int w, int h) :_w(w), _h(h), _pixels(h * w), _randoms(h * w) {
			XoroshiroPlus128 random;
			for (int i = 0; i < _randoms.size(); ++i) {
				_randoms[i] = random;
				random.jump();
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
		};
		const Pixel *pixel(int x, int y) const {
			return _pixels.data() + y * _w + x;
		}
		Pixel *pixel(int x, int y) {
			return _pixels.data() + y * _w + x;
		}

		PeseudoRandom *random(int x, int y) {
			return _randoms.data() + y * _w + x;
		}
	private:
		int _w = 0;
		int _h = 0;
		std::vector<Pixel> _pixels;
		std::vector<XoroshiroPlus128> _randoms;
	};

	inline double GTerm(const glm::dvec3 &p, double cosThetaP, const glm::dvec3 &q, double cosThetaQ) {
		return cosThetaP * cosThetaQ / glm::distance2(p, q);
	}

	inline bool has_value(const glm::dvec3 &c, double eps) {
		return glm::any(glm::greaterThanEqual(c, glm::dvec3(eps)));
	}

	///*
	//i = 0, p(0) = 1
	//i = 1, p(1) = exp(-1/k)
	//i = 2, p(2) = exp(-2/k)
	//*/
	//inline double russian_roulette_p(int i, double k = 40.0) {
	//	return std::exp(-double(i) / k);
	//}

	///*
	//i = 0, p(0) = 1
	//i = 1, p(0) * p(1) = p(-1/k)
	//i = 2, p(0) * p(1) * p(2) = p(-1/p) * p(-2/k) = exp(-3/k)
	//*/
	//inline double russian_roulette_p_cumulative(int i, double k = 40.0) {
	//	int sum = i * (1 + i) / 2;
	//	return std::exp(-double(sum) / k);
	//}

#define ENABLE_NEE 1
#define ENABLE_NEE_MIS 1
	inline glm::dvec3 radiance(const rt::SceneInterface &scene, glm::dvec3 ro, glm::dvec3 rd, PeseudoRandom *random) {
		// const double kSceneEPS = scene.adaptiveEps();
		const double kSceneEPS = 1.0e-6;
		const double kValueEPS = 1.0e-6;

		glm::dvec3 Lo;
		glm::dvec3 T(1.0);
		Material previous_m;
		double previous_pdf = 0.0;

		bool inside = false;

		constexpr int kDepth = 30;
		for (int i = 0; i < kDepth; ++i) {
			Material m;
			float tmin = 0.0f;
			glm::dvec3 wo = -rd;

			if (scene.intersect(ro, rd, &m, &tmin)) {
#if ENABLE_NEE
				// NEE
				// 最後のNEEは、PTと経路長をあわせるために、やらない
				if(i != (kDepth - 1)) {
					glm::dvec3 p = ro + rd * (double)tmin;
					for (auto it = scene.sampler_begin(); it != scene.sampler_end(); ++it) {
						glm::dvec3 q;
						glm::dvec3 n;
						glm::dvec3 Le;
						double pdf_area = 0.0;
						if ((*it)->can_sample(p) == false) {
							continue;
						}
						if (m->can_direct_sampling() == false) {
							continue;
						}

						(*it)->sample(random, p, &q, &n, &Le, &pdf_area);

						// 簡単なテスト
						// if (std::abs(pdf_area - (*it)->pdf_area(p, q)) > 1.0e-6) { abort(); }

						glm::dvec3 wi = glm::normalize(q - p);
						glm::dvec3 bxdf = m->bxdf(wo, wi);
						double cosThetaP = glm::dot(m->Ng, wi);
						double cosThetaQ = glm::dot(n, -wi);

						double g = GTerm(p, cosThetaP, q, cosThetaQ);

						glm::dvec3 contribution = T * bxdf * Le * g / pdf_area;

						if (has_value(contribution, kValueEPS)) {
							if (scene.occluded(p + m->Ng * kSceneEPS, q + n * kSceneEPS) == false) {
#if ENABLE_NEE_MIS
								double this_pdf = pdf_area;
								double other_pdf = m->pdf(wo, wi) * glm::dot(-n, wi) / glm::distance2(p, q);
								// double mis_weight = this_pdf / (this_pdf + other_pdf);
								double mis_weight = this_pdf * this_pdf / (this_pdf * this_pdf + other_pdf * other_pdf);
								Lo += contribution * mis_weight;
#else
								Lo += contribution;
#endif
							}
						}
					}
				}
#endif
				glm::dvec3 wi = m->sample(random, wo);
				glm::dvec3 bxdf = m->bxdf(wo, wi);
				glm::dvec3 emission = m->emission(wo);
				double pdf = m->pdf(wo, wi);
				double NoI = glm::dot(m->Ng, wi);
				double cosTheta = std::abs(NoI);

				if (inside) {
					T *= m->beers_law(tmin);
				}
				bool over_boundary = NoI < 0.0;
				if (over_boundary) {
					inside = !inside;
				}

				glm::dvec3 contribution = emission * T;

#if ENABLE_NEE_MIS
				if (has_value(contribution, kValueEPS)) {
					bool mis = false;
					if (i != 0) {
						if (auto sampler = m->direct_sampler()) {
							if (sampler->can_sample(previous_m->p) && previous_m->can_direct_sampling()) {
								double r = (double)tmin;
								double this_pdf = previous_pdf * glm::dot(m->Ng, wo) / (r * r);
								double other_pdf = sampler->pdf_area(previous_m->p, m->p);
								// double mis_weight = this_pdf * this_pdf / (this_pdf + other_pdf);
								double mis_weight = this_pdf * this_pdf / (this_pdf * this_pdf + other_pdf * other_pdf);
								mis = true;
								Lo += contribution * mis_weight;
							}
						}
					}

					if (mis == false) {
						Lo += contribution;
					}
				}
#elif ENABLE_NEE
				if (i == 0) {
					Lo += contribution;
				}
#else
				Lo += contribution;
#endif
				if (has_value(bxdf, kValueEPS)) {
					T *= bxdf * cosTheta / pdf;
				}
				else {
					break;
				}
				if (has_value(T, 1.0e-6) == false) {
					break;
				}

				// バイアスする方向は潜り込むときは逆転する
				ro = (ro + rd * (double)tmin) + (0.0 < NoI ? m->Ng : -m->Ng) * kSceneEPS;
				rd = wi;

				previous_pdf = pdf;
			}
			else {
				break;
			}
			previous_m = m;
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
			_badSampleCount = 0;
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
						PeseudoRandom *random = _image.random(x, y);
						glm::dvec3 o;
						glm::dvec3 d;
						_scene->camera.sampleRay(random, x, y, &o, &d);

						auto r = radiance(*_sceneInterface, o, d, random);

						bool badSample = false;
						for (int i = 0; i < r.length(); ++i) {
							if (glm::isfinite(r[i]) == false) {
								r[i] = 0.0;
								badSample = true;
							}
							if (r[i] < 0.0 || 1000.0 < r[i]) {
								r[i] = 0.0;
								badSample = true;
							}
						}
						_image.add(x, y, r);
						if (badSample) {
							_badSampleCount++;
						}
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

		int badSampleCount() const {
			return _badSampleCount.load();
		}
		std::shared_ptr<rt::Scene> _scene;
		std::shared_ptr<rt::SceneInterface> _sceneInterface;
		Image _image;
		int _steps = 0;
		std::atomic<int> _badSampleCount;
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
inline ofFloatPixels toOfLinear(const rt::Image &image) {
	ofFloatPixels pixels;
	pixels.allocate(image.width(), image.height(), OF_IMAGE_COLOR);
	float *dst = pixels.getPixels();

	for (int y = 0; y < image.height(); ++y) {
		for (int x = 0; x < image.width(); ++x) {
			int index = y * image.width() + x;
			const auto &px = *image.pixel(x, y);
			auto L = px.color / (double)px.sample;
			dst[index * 3 + 0] = L[0];
			dst[index * 3 + 1] = L[1];
			dst[index * 3 + 2] = L[2];
		}
	}
	return pixels;
}

std::shared_ptr<rt::Scene> scene;
std::shared_ptr<rt::PTRenderer> renderer;

bool isPowerOfTwo(uint32_t value)
{
	return value && !(value & (value - 1));
}

//--------------------------------------------------------------
void ofApp::setup() {
	_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
	_MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);

	ofxImGuiLite::initialize();

	ofSetVerticalSync(false);

	_camera.setNearClip(0.1);
	_camera.setFarClip(100.0);
	_camera.setDistance(5.0);

	rt::Stopwatch sw;

	scene = std::shared_ptr<rt::Scene>(new rt::Scene());
	rt::loadFromABC(ofToDataPath("cornelbox.abc").c_str(), *scene);
	// rt::loadFromABC(ofToDataPath("mitsuba.abc").c_str(), *scene);

	printf("load scene %f seconds\n", sw.elapsed());
	
	renderer = std::shared_ptr<rt::PTRenderer>(new rt::PTRenderer(scene));

	rt::CoupledBRDFConductor::load(
		ofToDataPath("baked/albedo_specular_conductor.bin").c_str(),
		ofToDataPath("baked/albedo_specular_conductor_avg.bin").c_str());
	rt::CoupledBRDFDielectrics::load(
		ofToDataPath("baked/albedo_specular_dielectrics.bin").c_str(),
		ofToDataPath("baked/albedo_specular_dielectrics_avg.bin").c_str());
}
void ofApp::exit() {
	ofxImGuiLite::shutdown();
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

	if (_showWireframe) {
		for (int i = 0; i < scene->geometries.size(); ++i) {
			static ofMesh mesh;
			mesh.clear();

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
			float tmin = 0.0f;
			if (renderer->sceneInterface().intersect(o, d, &m, &tmin)) {
				ofSetColor(255, 0, 0);
				auto p = o + d * (double)tmin;
				ofDrawLine(o.x, o.y, o.z, p.x, p.y, p.z);

				auto pn = p + m->Ng * 0.1;
				ofDrawLine(p.x, p.y, p.z, pn.x, pn.y, pn.z);
			}
			else {
				ofSetColor(255);
				auto p = o + d * 10.0;
				ofDrawLine(o.x, o.y, o.z, p.x, p.y, p.z);
			}
		}
	}

	if(_render) {
		renderer->step();

		ofDisableArbTex();

		if (ofGetFrameNum() % 5 == 0) {
			_image.setFromPixels(toOf(renderer->_image));
		}
		uint32_t n = renderer->stepCount();

		if (128 <= n && isPowerOfTwo(n)) {
			_image.setFromPixels(toOf(renderer->_image));
			char name[64];
			sprintf(name, "%dspp.png", n);
			_image.save(name);
			printf("elapsed %fs\n", ofGetElapsedTimef());
		}

		ofEnableArbTex();
	}

	_camera.end();

	ofDisableDepthTest();
	ofSetColor(255);

	ofxImGuiLite::ScopedImGui imgui;

	// camera control                                          for control clicked problem
	if (ImGui::IsWindowHovered(ImGuiHoveredFlags_AnyWindow) || (ImGui::IsAnyWindowFocused() && ImGui::IsAnyMouseDown())) {
		_camera.disableMouseInput();
	}
	else {
		_camera.enableMouseInput();
	}

	ImGui::SetNextWindowPos(ImVec2(20, 20), ImGuiSetCond_Appearing);
	ImGui::SetNextWindowSize(ImVec2(700, 600), ImGuiSetCond_Appearing);
	ImGui::SetNextWindowCollapsed(false, ImGuiSetCond_Appearing);
	ImGui::SetNextWindowBgAlpha(0.5f);

	ImGui::Begin("settings", nullptr);
	ImGui::Checkbox("render", &_render);
	ImGui::Checkbox("show wireframe", &_showWireframe);
	
	ImGui::Text("%d sample, fps = %.3f", renderer->stepCount(), ofGetFrameRate());
	ImGui::Text("%d bad sample", renderer->badSampleCount());
	if (_image.isAllocated()) {
		ofxImGuiLite::image(_image);
	}
	ImGui::End();
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

	if (key == 's') {
		ofFloatImage image = toOfLinear(renderer->_image);
		image.save("pt.exr");
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
