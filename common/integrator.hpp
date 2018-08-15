#pragma once

#include <atomic>
#include <tbb/tbb.h>
#include "scene_interface.hpp"

#define DEBUG_MODE 0

#define ENABLE_NEE 1
#define ENABLE_NEE_MIS 1

namespace rt {
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
	inline double GTerm(double cosThetaP, double cosThetaQ, double r2) {
		return cosThetaP * cosThetaQ / r2;
	}
	inline double GTerm(const glm::dvec3 &p, double cosThetaP, const glm::dvec3 &q, double cosThetaQ) {
		return GTerm(cosThetaP, cosThetaQ, glm::distance2(p, q));
	}

	inline bool has_value(const glm::dvec3 &c, double eps) {
		return glm::any(glm::greaterThanEqual(c, glm::dvec3(eps)));
	}

	inline glm::dvec3 radiance(const rt::SceneInterface &scene, glm::dvec3 ro, glm::dvec3 rd, PeseudoRandom *random) {
		const double kSceneEPS = scene.adaptiveEps();
		// const double kSceneEPS = 1.0e-6;
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
				if (i != (kDepth - 1)) {
					glm::dvec3 p = m->p;
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

						double pqDistance2 = glm::distance2(p, q);
						glm::dvec3 wi = (q - p) / std::sqrt(pqDistance2);
						
						double cosThetaP = glm::dot(m->Ng, wi);

						// 裏側に光源があるので早期棄却
						if (cosThetaP < 0.0) {
							continue;
						}

						// これはcan_sampleにおいてすでに裏面でないことが保証されている
						double cosThetaQ = glm::dot(n, -wi);

						glm::dvec3 bxdf = m->bxdf(wo, wi);

						double g = GTerm(cosThetaP, cosThetaQ, pqDistance2);

						glm::dvec3 contribution = T * bxdf * Le * g / pdf_area;

						if (has_value(contribution, kValueEPS)) {
							if (scene.occluded(p + m->Ng * kSceneEPS, q + n * kSceneEPS) == false) {
#if ENABLE_NEE_MIS
								double this_pdf = pdf_area;
								double other_pdf = m->pdf(wo, wi) * glm::dot(-n, wi) / pqDistance2;
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
				ro = m->p + (0.0 < NoI ? m->Ng : -m->Ng) * kSceneEPS;
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

	class PTRenderer {
	public:
		PTRenderer(std::shared_ptr<rt::Scene> scene)
			: _scene(scene)
			, _sceneInterface(new rt::SceneInterface(scene))
			, _image(scene->camera.imageWidth(), scene->camera.imageHeight()) {
			_badSampleNanCount = 0;
			_badSampleInfCount = 0;
			_badSampleNegativeCount = 0;
			_badSampleFireflyCount = 0;
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
					PeseudoRandom *random = _image.random(x, y);

					glm::dvec3 o;
					glm::dvec3 d;
					_scene->camera.sampleRay(random, x, y, &o, &d);

					auto r = radiance(*_sceneInterface, o, d, random);
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
							if (glm::isnan(r[i])) {
								_badSampleNanCount++;
								r[i] = 0.0;
							}
							else if (glm::isfinite(r[i]) == false) {
								_badSampleInfCount++;
								r[i] = 0.0;
							}
							else if (r[i] < 0.0) {
								_badSampleNegativeCount++;
								r[i] = 0.0;
							}
							if (10000.0 < r[i]) {
								_badSampleFireflyCount++;
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

		int badSampleNanCount() const {
			return _badSampleNanCount.load();
		}
		int badSampleInfCount() const {
			return _badSampleInfCount.load();
		}
		int badSampleNegativeCount() const {
			return _badSampleNegativeCount.load();
		}
		int badSampleFireflyCount() const {
			return _badSampleFireflyCount.load();
		}

		std::shared_ptr<rt::Scene> _scene;
		std::shared_ptr<rt::SceneInterface> _sceneInterface;
		Image _image;
		int _steps = 0;
		std::atomic<int> _badSampleNanCount;
		std::atomic<int> _badSampleInfCount;
		std::atomic<int> _badSampleNegativeCount;
		std::atomic<int> _badSampleFireflyCount;
	};
}