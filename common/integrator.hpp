#pragma once

#include <atomic>
#include <tbb/tbb.h>
#include "scene_interface.hpp"
#include "online.hpp"

#define DEBUG_MODE 0
#define ENABLE_ADAPTIVE_SAMPLING 0

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
			if (_pixels[index].sample & 0x1) {
				_pixels[index].color_half += c;
			}
		}

		struct Pixel {
			int sample = 0;
			glm::dvec3 color;
			glm::dvec3 color_half;

			double ep() const {
				if (sample & 0x1) {
					abort();
				}
				glm::dvec3 I = color / double(sample);
				glm::dvec3 A = 2.0 * color_half / double(sample);
				glm::dvec3 d = glm::abs(I - A);
				double denom = std::sqrt(I.x + I.y + I.z);
				if (denom < 1.0e-12) {
					return 0.0;
				}
				double e = (d.x + d.y + d.z) / denom;
				return e;
			}
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
				auto nee = [&]() {
					if (i == (kDepth - 1)) {
						return;
					}
					if (m->can_direct_sampling() == false) {
						return;
					}
					glm::dvec3 p = m->p;
					
					static thread_local std::vector<double> importance_storage;
					LightSelector selector(p, scene.sampler_begin(), scene.sampler_end(), importance_storage);
					if (selector.can_sample() == false) {
						return;
					}

					double p_choice = 0.0;
					auto sampler = selector.choice(random, &p_choice);
					glm::dvec3 q;
					glm::dvec3 n;
					glm::dvec3 Le;
					double pdf_area = 0.0;
					sampler->sample(random, p, &q, &n, &Le, &pdf_area);

					double pqDistance2 = glm::distance2(p, q);
					glm::dvec3 wi = (q - p) / std::sqrt(pqDistance2);

					double cosThetaP = glm::dot(m->Ng, wi);

					// 裏側に光源があるので早期棄却
					if (cosThetaP < 0.0) {
						return;
					}

					// これはcan_sampleにおいてすでに裏面でないことが保証されている
					double cosThetaQ = glm::dot(n, -wi);

					glm::dvec3 bxdf = m->bxdf(wo, wi);

					double g = GTerm(cosThetaP, cosThetaQ, pqDistance2);

					glm::dvec3 contribution = T * bxdf * Le * g / pdf_area / p_choice;

					if (has_value(contribution, kValueEPS) == false) {
						return;
					}
					if (scene.occluded(p + m->Ng * kSceneEPS, q + n * kSceneEPS)) {
						return;
					}
#if ENABLE_NEE_MIS
					double this_pdf = pdf_area * p_choice;
					double other_pdf = m->pdf(wo, wi) * glm::dot(-n, wi) / pqDistance2;
					// double mis_weight = this_pdf / (this_pdf + other_pdf);
					double mis_weight = this_pdf * this_pdf / (this_pdf * this_pdf + other_pdf * other_pdf);
					Lo += contribution * mis_weight;
#else
					Lo += contribution;
#endif
				};
				nee();
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
								static thread_local std::vector<double> importance_storage;
								LightSelector selector(previous_m->p, scene.sampler_begin(), scene.sampler_end(), importance_storage);

								double r = (double)tmin;
								double this_pdf = previous_pdf * glm::dot(m->Ng, wo) / (r * r);
								double other_pdf = sampler->pdf_area(previous_m->p, m->p) * selector.p(sampler);
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

	// x0 <= u < x1
	// y0 <= v < y1
	class Region {
	public:
		int x0 = 0;
		int y0 = 0;
		int x1 = 0;
		int y1 = 0;

		int area() const {
			return (x1 - x0) * (y1 - y0);
		}
		std::tuple<Region, Region> divide() const {
			int w = x1 - x0;
			int h = y1 - y0;
			if (w < h) {
				int cy = y0 + h / 2;
				Region L;
				L.x0 = x0;
				L.y0 = y0;
				L.x1 = x1;
				L.y1 = cy;

				Region R;
				R.x0 = x0;
				R.y0 = cy;
				R.x1 = x1;
				R.y1 = y1;
				return std::tuple<Region, Region>(L, R);
			}
			else {
				int cx = x0 + w / 2;
				Region L;
				L.x0 = x0;
				L.y0 = y0;
				L.x1 = cx;
				L.y1 = y1;

				Region R;
				R.x0 = cx;
				R.y0 = y0;
				R.x1 = x1;
				R.y1 = y1;
				return std::tuple<Region, Region>(L, R);
			}
			return std::tuple<Region, Region>();
		}
		void divide(int min_area, std::vector<Region> &regions) const {
			rt::Region A, B;
			std::tie(A, B) = divide();
			if (min_area < A.area() && min_area < B.area()) {
				A.divide(min_area, regions);
				B.divide(min_area, regions);
			}
			else {
				regions.push_back(*this);
			}
		}
	};

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

			std::vector<rt::Region> regions;
			rt::Region whole;
			whole.x0 = 0;
			whole.y0 = 0;
			whole.x1 = _image.width();
			whole.y1 = _image.height();
			whole.divide(32 * 32, regions);

			for (int i = 0; i < regions.size(); ++i) {
				Block b;
				b.region = regions[i];
				_blocks.push_back(b);
			}
		}
		void step() {
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
			
#if ENABLE_ADAPTIVE_SAMPLING
			auto step_block = [&](const Block &block){
				const Region &region = block.region;
				for (int y = region.y0; y < region.y1; ++y) {
					for (int x = region.x0; x < region.x1; ++x) {
						PeseudoRandom *random = _image.random(x, y);
						glm::dvec3 o;
						glm::dvec3 d;
						_scene->camera.sampleRay(random, x, y, &o, &d);

						auto r = radiance(*_sceneInterface, o, d, random);

						check_sample(r);

						_image.add(x, y, r);
					}
				}
			};
			auto update_metric = [&](Block &block) {
				const Region &region = block.region;
				rt::OnlineMean<double> ed;
				for (int y = region.y0; y < region.y1; ++y) {
					for (int x = region.x0; x < region.x1; ++x) {
						int index = y * _image.width() + x;
						const auto &px = *_image.pixel(x, y);
						ed.addSample(px.ep());
					}
				}
				block.metric = ed.mean();
			};

			const int kAdaptiveBegin = 32;
			const int kAdaptiveFreq = 8;
			const int kActiveBlocks = _blocks.size() / 3;

			if (_steps < kAdaptiveBegin) {
				tbb::parallel_for(tbb::blocked_range<int>(0, _blocks.size()), [&](const tbb::blocked_range<int> &range) {
					for (int i = range.begin(); i < range.end(); ++i) {
						step_block(_blocks[i]);
						_blocks[i].sample++;
					}
				});
			}
			else {
				// kAdaptiveBegin <= _steps
				if ((_steps - kAdaptiveBegin) % kAdaptiveFreq == 0) {
					int metric_blocks = _steps == kAdaptiveBegin ? _blocks.size() : kActiveBlocks;
					tbb::parallel_for(tbb::blocked_range<int>(0, metric_blocks), [&](const tbb::blocked_range<int> &range) {
						for (int i = range.begin(); i < range.end(); ++i) {
							update_metric(_blocks[i]);
						}
					});
					std::sort(_blocks.begin(), _blocks.end(), [](const Block &a, const Block &b) {
						return a.metric > b.metric;
					});
				}
				// kActiveBlocks
				// _blocks.size()
				tbb::parallel_for(tbb::blocked_range<int>(0, kActiveBlocks), [&](const tbb::blocked_range<int> &range) {
					for (int i = range.begin(); i < range.end(); ++i) {
						step_block(_blocks[i]);
						_blocks[i].sample++;
					}
				});
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

						check_sample(r);

						_image.add(x, y, r);
					}
				}
			});
#endif

#endif
			_steps++;
		}
		int stepCount() const {
			return _steps;
		}

		const rt::SceneInterface &sceneInterface() const {
			return *_sceneInterface;
		}

		void check_sample(glm::dvec3 &r) {
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

		struct Block {
			Region region;
			double metric = 0.0;
			int sample = 0;
		};
		std::vector<Block> _blocks;
	};
}