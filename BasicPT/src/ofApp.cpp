#include "alembic_loader.hpp"
#include "render_object.hpp"
#include "microfacet.hpp"
#include "bicubic.hpp"
#include "material.hpp"
#include "value_prportional_sampler.hpp"
#include "geometry.hpp"
#include "stopwatch.hpp"

#include <random>
#include <xmmintrin.h>
#include <pmmintrin.h>
#include <tbb/tbb.h>
#include <embree3/rtcore.h>

#include "ofApp.h"

#define DEBUG_MODE 0

namespace rt {
	struct Rectangle {
		Rectangle() {}
		Rectangle(glm::dvec3 s, glm::dvec3 ex, glm::dvec3 ey)
		:_s(s)
		,_ex(ex)
		,_ey(ey) {
			_exLength = glm::length(ex);
			_eyLength = glm::length(ey);
			_x = ex / _exLength;
			_y = ey / _exLength;
			_z = glm::cross(_x, _y);
		}
		glm::dvec3 sample(double u, double v) const {
			return _s + _ex * u + _ey * v;
		}
		double area() const {
			return _exLength * _eyLength;
		}
		glm::dvec3 normal() const {
			return _z;
		}
		glm::dvec3 s() const {
			return _s;
		}
		glm::dvec3 ex() const {
			return _ex;
		}
		glm::dvec3 ey() const {
			return _ey;
		}
		glm::dvec3 x() const {
			return _x;
		}
		glm::dvec3 y() const {
			return _y;
		}
		glm::dvec3 z() const {
			return _z;
		}
		double exLength() const {
			return _exLength;
		}
		double eyLength() const {
			return _eyLength;
		}
		glm::dvec3 _s;
		glm::dvec3 _ex;
		glm::dvec3 _ey;
		double _exLength = 0.0;
		double _eyLength = 0.0;
		glm::dvec3 _x;
		glm::dvec3 _y;
		glm::dvec3 _z;
	};

	class SphericalRectangleSamplerCoordinate {
	public:
		SphericalRectangleSamplerCoordinate(const Rectangle &rectangle, const glm::dvec3 &o):_rectangle(rectangle) {
			_o = o;

			glm::dvec3 d = rectangle.s() - o;
			_x0 = glm::dot(d, rectangle.x());
			_y0 = glm::dot(d, rectangle.y());
			_z0 = glm::dot(d, rectangle.z());

			_x1 = _x0 + rectangle.exLength();
			_y1 = _y0 + rectangle.eyLength();

			// z flip
			if (_z0 > 0.0) {
				_z0 *= -1.0;
				_rectangle._z *= -1.0;
			}

			_v[0][0] = glm::vec3(_x0, _y0, _z0);
			_v[0][1] = glm::vec3(_x0, _y1, _z0);
			_v[1][0] = glm::vec3(_x1, _y0, _z0);
			_v[1][1] = glm::vec3(_x1, _y1, _z0);

			auto calculate_n = [](glm::vec3 a, glm::vec3 b) {
				glm::vec3 c = glm::cross(a, b);
				return glm::normalize(c);
			};

			_n[0] = calculate_n(_v[0][0], _v[1][0]);
			_n[1] = calculate_n(_v[1][0], _v[1][1]);
			_n[2] = calculate_n(_v[1][1], _v[0][1]);
			_n[3] = calculate_n(_v[0][1], _v[0][0]);

			_g[0] = std::acos(glm::dot(-_n[0], _n[1]));
			_g[1] = std::acos(glm::dot(-_n[1], _n[2]));
			_g[2] = std::acos(glm::dot(-_n[2], _n[3]));
			_g[3] = std::acos(glm::dot(-_n[3], _n[0]));
		}

		double solidAngle() const {
			return _g[0] + _g[1] + _g[2] + _g[3] - glm::two_pi<double>();
		}

		double minPhiU() const {
			return -_g[2] - _g[3] + glm::two_pi<double>();
		}
		double maxPhiU() const {
			return solidAngle() - _g[2] - _g[3] + glm::two_pi<double>();
		}
		glm::dvec3 sample(double u, double v) const {
			double AQ = solidAngle();

			double phi_u = u * AQ - _g[2] - _g[3] + glm::two_pi<double>();

			auto safeSqrt = [](double x) {
				return std::sqrt(std::max(x, 0.0));
			};

			double b0 = _n[0].z;
			double b1 = _n[2].z;
			double fu = (std::cos(phi_u) * b0 - b1) / std::sin(phi_u);
			double cu = std::copysign(1.0, fu) / safeSqrt(fu * fu + b0 * b0);
			double xu = -cu * _z0 / safeSqrt(1.0 - cu * cu);

			double d = std::sqrt(xu * xu + _z0 * _z0);
			double h0 = _y0 / std::sqrt(d * d + _y0 * _y0);
			double h1 = _y1 / std::sqrt(d * d + _y1 * _y1);
			double hv = glm::mix(h0, h1, v);
			double yv = hv * d / safeSqrt(1.0 - hv * hv);
			return _o + xu * _rectangle.x() + yv * _rectangle.y() + _z0 * _rectangle.z();
		}

		Rectangle _rectangle;

		glm::dvec3 _o;

		// local coord : 
		double _x0;
		double _x1;
		double _y0;
		double _y1;
		double _z0;

		glm::dvec3 _v[2][2];
		glm::dvec3 _n[4];
		double _g[4];
	};
}

namespace rt {
	/*
	片面の場合、サンプリングは迷わないが、法線により衝突判定は早期に棄却可能
	両面ある場合、明示的なサンプリングする時、裏面をサンプリングするのは積分範囲の無駄であるので、
	見えている方向だけをサンプリングできる。
	*/
	class IDirectSampler {
	public:
		virtual double pdf_area(glm::dvec3 o, glm::dvec3 p) const = 0;
		virtual void sample(PeseudoRandom *random, glm::dvec3 o, glm::dvec3 *p, glm::dvec3 *n, glm::dvec3 *Le, double *pdf_area) const = 0;
	};

	class TriangleAreaSampler : public IDirectSampler {
	public:
		TriangleAreaSampler(glm::dvec3 a, glm::dvec3 b, glm::dvec3 c, bool doubleSided, glm::dvec3 Le)
		:_a(a)
		,_b(b)
		,_c(c)
		,_doubleSided(doubleSided)
		,_Le(Le) {
			_n = triangleNormal(_a, _b, _c);
			_one_over_area = 1.0 / triangleArea(_a, _b, _c);
		}

		virtual double pdf_area(glm::dvec3 o, glm::dvec3 p) const {
			// 裏面はサンプリングされないため、確率密度は0
			if (_doubleSided == false) {
				glm::dvec3 d = p - o;
				bool backfacing = 0.0 < glm::dot(d, _n);
				if (backfacing) {
					return 0.0;
				}
			}
			return _one_over_area;
		}
		
		virtual void sample(PeseudoRandom *random, glm::dvec3 o, glm::dvec3 *p, glm::dvec3 *n, glm::dvec3 *Le, double *pdf_area) const {
			*p = uniform_on_triangle(random->uniform(), random->uniform()).evaluate(_a, _b, _c);

			glm::dvec3 d = *p - o;
			bool backfacing = 0.0 < glm::dot(d, _n);

			if (_doubleSided) {
				*n = backfacing ? -_n : _n;
				*Le = _Le;
				*pdf_area = _one_over_area;
			}
			else {
				*n = _n;
				*Le = backfacing ? glm::dvec3(0.0) : _Le;

				// 裏面はサンプリングされないため、確率密度は0
				*pdf_area = backfacing ? 0.0 : _one_over_area;
			}
		}
	private:
		glm::dvec3 _Le;
		glm::dvec3 _a, _b, _c;
		glm::dvec3 _n;
		bool _doubleSided = false;
		double _one_over_area = 0.0;
	};
	class SphericalRectangleSampler : public IDirectSampler {
	public:
		SphericalRectangleSampler(glm::dvec3 s, glm::dvec3 ex, glm::dvec3 ey, bool doubleSided, glm::dvec3 Le)
		:_doubleSided(doubleSided)
		,_Le(Le)
		,_q(s, ex, ey)
		{

		}

		virtual void sample(PeseudoRandom *random, glm::dvec3 o, glm::dvec3 *p, glm::dvec3 *n, glm::dvec3 *Le, double *pdf_area) const
		{
			SphericalRectangleSamplerCoordinate sampler(_q, o);
			*p = sampler.sample(random->uniform(), random->uniform());
			glm::dvec3 d = *p - o;
			bool backfacing = 0.0 < glm::dot(d, _q.normal());

			if (_doubleSided) {
				*n = backfacing ? -_q.normal() : _q.normal();
				*Le = _Le;
			}
			else {
				*n = _q.normal();
				*Le = backfacing ? glm::dvec3(0.0) : _Le;
			}

			if (_doubleSided == false && backfacing) {
				*pdf_area = 0.0;
			} else {
				double dLength2 = glm::length2(d);
				double cosTheta = glm::dot(-*n, d / std::sqrt(dLength2));
				double pw = 1.0 / sampler.solidAngle();
				*pdf_area = pw * cosTheta / dLength2;
			}
		}
		virtual double pdf_area(glm::dvec3 o, glm::dvec3 p) const {
			glm::dvec3 d = p - o;
			bool backfacing = 0.0 < glm::dot(d, _q.normal());

			// 裏面はサンプリングされないため、確率密度は0
			if (_doubleSided == false) {
				if (backfacing) {
					return 0.0;
				}
			}
			SphericalRectangleSamplerCoordinate sampler(_q, o);

			glm::dvec3 n = backfacing ? -_q.normal() : _q.normal();
			double dLength2 = glm::length2(d);
			double cosTheta = glm::dot(-n, d / std::sqrt(dLength2));
			double pw = 1.0 / sampler.solidAngle();
			return pw * cosTheta / glm::distance2(o, p);
		}
	private:
		glm::dvec3 _Le;
		bool _doubleSided = false;
		Rectangle _q;
	};

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
				const Geometry& g = _scene->geometries[i];
				for (int j = 0; j < g.primitives.size(); ++j) {
					const Geometry::Primitive &p = g.primitives[j];
					const IMaterial *m = g.primitives[j].material.get();
					const LambertianMaterial *lambertian = dynamic_cast<const LambertianMaterial *>(m);
					if (lambertian && lambertian->isEmission()) {
						if (auto sample = lambertian->samplingStrategy.get<AreaSample>()) {
							auto a = g.points[p.indices[0]].P;
							auto b = g.points[p.indices[1]].P;
							auto c = g.points[p.indices[2]].P;
							TriangleAreaSampler *sampler = new TriangleAreaSampler(a, b, c, lambertian->backEmission, lambertian->Le);
							_directSamplers.emplace_back(sampler);
						}
						else if (auto sample = lambertian->samplingStrategy.get<SphericalRectangleSample>())
						{
							if (sample->triangleIndex == 0) {
								SphericalRectangleSampler *sampler = new SphericalRectangleSampler(sample->s, sample->ex, sample->ey, lambertian->backEmission, lambertian->Le);
								_directSamplers.emplace_back(sampler);
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
			const auto &prim = _scene->geometries[index].primitives[rayhit.hit.primID];
			*material = prim.material;

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
			return true;
		}

		const Camera &camera() const {
			return _scene->camera;
		}

		std::vector<std::shared_ptr<IDirectSampler>>::const_iterator sampler_begin() const {
			return _directSamplers.begin();
		}
		std::vector<std::shared_ptr<IDirectSampler>>::const_iterator sampler_end() const {
			return _directSamplers.end();
		}
		int samplerCount() const {
			return _directSamplers.size();
		}

		std::shared_ptr<rt::Scene> _scene;
		RTCDevice _embreeDevice = nullptr;
		RTCScene _embreeScene = nullptr;
		mutable RTCIntersectContext _context;

		std::vector<std::shared_ptr<IDirectSampler>> _directSamplers;

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

#define ENABLE_NEE 1
	inline glm::dvec3 radiance(const rt::SceneInterface &scene, glm::dvec3 ro, glm::dvec3 rd, PeseudoRandom *random) {
		const double kSceneEPS = scene.adaptiveEps();

		glm::dvec3 Lo;
		glm::dvec3 T(1.0);
		for (int i = 0; i < 10; ++i) {
			Material m;
			float tmin = 0.0f;
			glm::dvec3 wo = -rd;

			if (scene.intersect(ro, rd, &m, &tmin)) {
#if ENABLE_NEE
				// NEE
				{
					glm::dvec3 p = ro + rd * (double)tmin;

					// 1 sample only
					int samplerCount = scene.samplerCount();
					int sample_index = (int)random->uniform(0.0, samplerCount);
					sample_index = std::min(sample_index, samplerCount - 1);
					double sampler_select_p = 1.0 / samplerCount;

					for (auto it = scene.sampler_begin(); it != scene.sampler_end(); ++it) {
						if (std::distance(scene.sampler_begin(), it) != sample_index) {
							continue;
						}

						glm::dvec3 q;
						glm::dvec3 n;
						glm::dvec3 Le;
						double pdf_area = 0.0;
						(*it)->sample(random, p, &q, &n, &Le, &pdf_area);

						// 簡単なテスト
						// if (std::abs(pdf_area - (*it)->pdf_area(p, q)) > 1.0e-6) { abort(); }

						glm::dvec3 wi = glm::normalize(q - p);
						glm::dvec3 bxdf = m->bxdf(wo, wi);
						double cosThetaP = glm::abs(glm::dot(m->Ng, wi));
						double cosThetaQ = glm::dot(n, -wi);

						double g = GTerm(p, cosThetaP, q, cosThetaQ);

						glm::dvec3 contribution = T * bxdf * Le * g;

						/* 裏面は発光しない */
						// 0.0 < cosThetaQ && 
						if (glm::any(glm::greaterThanEqual(contribution, glm::dvec3(glm::epsilon<double>())))) {
							if (scene.occluded(p + m->Ng * kSceneEPS, q + n * kSceneEPS) == false) {
								Lo += contribution / (pdf_area * sampler_select_p);
							}
						}
					}
				}
#endif
				glm::dvec3 wi = m->sample(random, wo);
				glm::dvec3 bxdf = m->bxdf(wo, wi);
				glm::dvec3 emission = m->emission(wo);
				double pdf = m->pdf(wo, wi);
				double cosTheta = std::abs(glm::dot(m->Ng, wi));


#if ENABLE_NEE
				if (i == 0) {
					Lo += emission * T;
				}
#else
				Lo += emission * T;
#endif
				if (glm::any(glm::greaterThanEqual(bxdf, glm::dvec3(1.0e-6f)))) {
					T *= bxdf * cosTheta / pdf;
				}
				else {
					break;
				}

				ro = (ro + rd * (double)tmin) + m->Ng * kSceneEPS;
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
						PeseudoRandom *random = _image.random(x, y);
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

	if (0) {
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

	{
		renderer->step();

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
		_image.draw(30, 30);
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

	//if (key == 's') {
	//	_image.save("pt.png");
	//	ofFloatImage img(_image);
	//}
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
