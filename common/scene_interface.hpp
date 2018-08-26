#pragma once
#include <embree3/rtcore.h>

#include "render_object.hpp"
#include "microfacet.hpp"
#include "bicubic.hpp"
#include "material.hpp"
#include "value_prportional_sampler.hpp"
#include "geometry.hpp"
#include "stopwatch.hpp"
#include "direct_sampler.hpp"

namespace rt {
	inline void EmbreeErorrHandler(void* userPtr, RTCError code, const char* str) {
		printf("Embree Error [%d] %s\n", code, str);
	}
	class SceneInterface {
	public:
		SceneInterface(std::shared_ptr<rt::Scene> scene) :_scene(scene) {
			_embreeDevice = rtcNewDevice("set_affinity=1");
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

			glm::dvec3 Ng = prim.Ng;

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

		std::vector<IDirectSampler *>::const_iterator sampler_begin() const {
			return _directSamplers.begin();
		}
		std::vector<IDirectSampler *>::const_iterator sampler_end() const {
			return _directSamplers.end();
		}
		int samplerCount() const {
			return _directSamplers.size();
		}

		std::shared_ptr<rt::Scene> _scene;
		RTCDevice _embreeDevice = nullptr;
		RTCScene _embreeScene = nullptr;
		mutable RTCIntersectContext _context;

		std::vector<IDirectSampler *> _directSamplers;

		double _sceneAdaptiveEps = 0.0f;
	};
}