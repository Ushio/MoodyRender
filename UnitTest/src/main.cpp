#include "ofMain.h"

#define CATCH_CONFIG_RUNNER
#include <catch.hpp>

#include <set>
#include <memory>

#include "online.hpp"
#include "peseudo_random.hpp"
#include "composite_simpson.hpp"
#include "adaptive_simpson.hpp"
#include "simpson_helper.hpp"
#include "coordinate.hpp"
#include "microfacet.hpp"
#include "serializable_buffer.hpp"
#include "material.hpp"
#include "geometry.hpp"
#include "randomsampler.hpp"

TEST_CASE("online", "[online]") {
	SECTION("online") {
		rt::Xor64 random;
		for (int i = 0; i < 100; ++i) {
			std::vector<double> xs;
			rt::OnlineVariance<double> ov;
			for (int j = 0; j < 100; ++j) {
				double x = random.uniform(0.0, 10.0);
				xs.push_back(x);
				ov.addSample(x);

				double mean;
				double variance;
				rt::mean_and_variance(xs, &mean, &variance);

				REQUIRE(std::abs(mean - ov.mean()) < 1.0e-9);
				REQUIRE(std::abs(variance - ov.variance()) < 1.0e-9);
			}
		}
	}
}

TEST_CASE("simpson", "[simpson]") {
	// example is from:
	//     http://mathfaculty.fullerton.edu/mathews/n2003/AdaptiveQuadMod.html

	SECTION("duplicated sample detection adaptive_simpson") {
		auto f = [](double x) {
			return (14.0 * x - 11.0 * x * x) * std::exp(-2.0 * x);
		};
		std::set<double> evaluated;
		rt::SimpsonRange<double> simpson([&](double x) {
			double v = f(x);
			REQUIRE(evaluated.count(v) == 0);
			evaluated.insert(v);
			return v;
		}, 0.0, 10.0);
		rt::adaptive_simpson(simpson, 1.0e-13);
	}

	SECTION("example A") {
		// http://www.wolframalpha.com/input/?i=Plot%5B(14x+-+11x%5E2)+exp(-2x),+%7Bx,+0,+5%7D%5D
		auto f = [](double x) {
			return (14.0 * x - 11.0 * x * x) * std::exp(-2.0 * x);
		};
		auto f_integral_analytic = [](double a, double b) {
			auto f = [](double x) {
				// http://www.wolframalpha.com/input/?i=Integral%5B(14x+-+11x%5E2)+exp(-2x),+x%5D
				return 0.25 * std::exp(-2.0 * x) * (22.0 * x * x - 6.0 * x - 3.0);
			};
			return f(b) - f(a);
		};

		rt::Xor64 random;

		SECTION("composite_simpson") {
			for (int i = 0; i < 10000; ++i) {
				double a = random.uniform(0.0, 5.0);
				double b = random.uniform(0.0, 5.0);
				double numerical = rt::composite_simpson<double>(f, 300, a, b);
				double analytical = f_integral_analytic(a, b);
				double diff = std::abs(numerical - analytical);
				CAPTURE(a);
				CAPTURE(b);
				CAPTURE(diff);
				REQUIRE(diff < 1.0e-6);
			}
		}

		SECTION("adaptive_simpson") {
			// rt::OnlineVariance<double> ov;
			for (int i = 0; i < 10000; ++i) {
				double a = random.uniform(0, 5.0);
				double b = random.uniform(0, 5.0);

				int sample = 0;
				rt::SimpsonRange<double> simpson([&](double x) {
					sample++;
					return f(x);
				}, a, b);

				double numerical = rt::adaptive_simpson(simpson, 1.0e-10);
				double analytical = f_integral_analytic(a, b);
				double diff = std::abs(numerical - analytical);
				CAPTURE(a);
				CAPTURE(b);
				CAPTURE(sample);
				REQUIRE(diff < 1.0e-6);

				// ov.addSample(sample);
			}
			// printf("sample mean: %f, sd: %f\n", ov.mean(), std::sqrt(ov.variance()));
		}
	}
}

TEST_CASE("LambertianMaterial", "[LambertianMaterial]") {

	SECTION("white furnance LambertianMaterial") {
		using namespace rt;

		rt::Xor64 *random = new rt::Xor64();
		for (int j = 0; j < 100; ++j) {

			// alphaが小さい場合、simpsonによる積分が適さない
			double alpha = random->uniform(0.1, 1.0);
			glm::dvec3 Ng(0.0, 0.0, 1.0);
			glm::dvec3 wo = LambertianSampler::sample(random, Ng);

			LambertianMaterial m;
			m.Ng = Ng;
			m.Le = glm::dvec3(0.0);
			m.R = glm::dvec3(1.0);

			double result = hemisphere_composite_simpson<double>([&](double theta, double phi) {
				glm::dvec3 wi = rt::polar_to_cartesian((double)theta, (double)phi);

				glm::dvec3 brdf = m.bxdf(wo, wi);
				double cosTheta = glm::dot(m.Ng, wi);

				REQUIRE(std::abs(brdf.x - brdf.y) < 1.0e-6);
				REQUIRE(std::abs(brdf.y - brdf.z) < 1.0e-6);

				return brdf.r * cosTheta;
			}, 500);

			CAPTURE(result);
			CAPTURE(alpha);
			CAPTURE(glm::dot(Ng, wo));
			REQUIRE(std::abs(result - 1.0) < 1.0e-2);
		}
	}
}

TEST_CASE("microfacet", "[microfacet]") {
	rt::CoupledBRDFConductor::load(
		ofToDataPath("baked/albedo_specular_conductor.bin").c_str(),
		ofToDataPath("baked/albedo_specular_conductor_avg.bin").c_str());
	rt::CoupledBRDFDielectrics::load(
		ofToDataPath("baked/albedo_specular_dielectrics.bin").c_str(),
		ofToDataPath("baked/albedo_specular_dielectrics_avg.bin").c_str());

	SECTION("hemisphere_composite_simpson") {
		double result = rt::hemisphere_composite_simpson<double>([](double theta, double phi) {
			return 1.0 / glm::pi<double>() * std::cos(theta);
		}, 100);
		REQUIRE(std::abs(result - 1.0) < 1.0e-8);
	}

	SECTION("beckmann normalization") {
		using namespace rt;

		rt::Xor64 random;
		for (int j = 0; j < 32; ++j) {
			// alphaが小さい場合、simpsonによる積分が適さない
			double alpha = random.uniform(0.1, 1.0);

			double result = hemisphere_composite_simpson<double>([&](double theta, double phi) {
				glm::dvec3 wi = rt::polar_to_cartesian((double)theta, (double)phi);
				glm::dvec3 Ng(0.0, 0.0, 1.0);
				glm::dvec3 sample_m = rt::polar_to_cartesian((double)theta, (double)phi);
				double cosTheta = glm::dot(Ng, sample_m);
				double value = rt::D_Beckmann(Ng, sample_m, alpha) * cosTheta;
				return (double)value;
			}, 500);
			CAPTURE(result);
			CAPTURE(alpha);
			REQUIRE(std::abs(result - 1.0) < 1.0e-4);

		}
	}

	SECTION("white furnance test MicrofacetCoupledConductorMaterial") {
		using namespace rt;

		rt::Xor64 *random = new rt::Xor64();
		for (int j = 0; j < 32; ++j) {

			// alphaが小さい場合、simpsonによる積分が適さない
			double alpha = random->uniform(0.1, 1.0);
			glm::dvec3 Ng(0.0, 0.0, 1.0);
			glm::dvec3 wo = LambertianSampler::sample(random, Ng);

			MicrofacetCoupledConductorMaterial m;
			m.Ng = Ng;
			m.alpha = alpha;
			m.useFresnel = false;

			double result = hemisphere_composite_simpson<double>([&](double theta, double phi) {
				glm::dvec3 wi = rt::polar_to_cartesian((double)theta, (double)phi);

				glm::dvec3 brdf = m.bxdf(wo, wi);
				double cosTheta = glm::dot(m.Ng, wi);

				REQUIRE(std::abs(brdf.x - brdf.y) < 1.0e-6);
				REQUIRE(std::abs(brdf.y - brdf.z) < 1.0e-6);

				return brdf.r * cosTheta;
			}, 500);

			CAPTURE(result);
			CAPTURE(alpha);
			CAPTURE(glm::dot(Ng, wo));
			REQUIRE(std::abs(result - 1.0) < 1.0e-2);

			printf("coupled conductor %.10f\n", result);
		}
	}

	SECTION("white furnance test MicrofacetCoupledDielectricsMaterial") {
		using namespace rt;

		rt::Xor64 *random = new rt::Xor64();
		for (int j = 0; j < 32; ++j) {

			// alphaが小さい場合、simpsonによる積分が適さない
			double alpha = random->uniform(0.1, 1.0);
			glm::dvec3 Ng(0.0, 0.0, 1.0);
			glm::dvec3 wo = LambertianSampler::sample(random, Ng);

			MicrofacetCoupledDielectricsMaterial m;
			m.Ng = Ng;
			m.alpha = alpha;

			double result = hemisphere_composite_simpson<double>([&](double theta, double phi) {
				glm::dvec3 wi = rt::polar_to_cartesian((double)theta, (double)phi);

				glm::dvec3 brdf = m.bxdf(wo, wi);
				double cosTheta = glm::dot(m.Ng, wi);

				REQUIRE(std::abs(brdf.x - brdf.y) < 1.0e-6);
				REQUIRE(std::abs(brdf.y - brdf.z) < 1.0e-6);

				return brdf.r * cosTheta;
			}, 500);

			CAPTURE(result);
			CAPTURE(alpha);
			CAPTURE(glm::dot(Ng, wo));
			REQUIRE(std::abs(result - 1.0) < 1.0e-2);
			printf("coupled dielectrics %.10f\n", result);
		}
	}

	SECTION("white furnance test MC MicrofacetCoupledConductorMaterial") {
		using namespace rt;

		rt::Xor64 *random = new rt::Xor64();
		for (int j = 0; j < 32; ++j) {

			// alphaが小さい場合、simpsonによる積分が適さない
			double alpha = random->uniform(0.1, 1.0);
			glm::dvec3 Ng(0.0, 0.0, 1.0);
			glm::dvec3 wo = LambertianSampler::sample(random, Ng);

			MicrofacetCoupledConductorMaterial m;
			m.Ng = Ng;
			m.alpha = alpha;
			m.useFresnel = false;

			OnlineMean<double> mean;

			for (int i = 0; i < 500000; ++i) {
				glm::dvec3 wi = m.sample(random, wo);
				glm::dvec3 bxdf = m.bxdf(wo, wi);
				double pdf = m.pdf(wo, wi);
				double cosTheta = glm::dot(m.Ng, wi);

				REQUIRE(std::abs(bxdf.x - bxdf.y) < 1.0e-6);
				REQUIRE(std::abs(bxdf.y - bxdf.z) < 1.0e-6);

				glm::dvec3 value;
				if (glm::any(glm::greaterThanEqual(bxdf, glm::dvec3(1.0e-6f)))) {
					value = bxdf * cosTheta / pdf;
				}
				else {
					value = glm::dvec3(0.0);
				}

				mean.addSample(value.x);
			}
			double result = mean.mean();

			CAPTURE(alpha);
			CAPTURE(glm::dot(Ng, wo));
			REQUIRE(std::abs(result - 1.0) < 1.0e-2);
			printf("MC coupled conductor %.10f\n", result);
		}
	}

	SECTION("white furnance test MC MicrofacetCoupledDielectricsMaterial") {
		using namespace rt;

		rt::Xor64 *random = new rt::Xor64();
		for (int j = 0; j < 32; ++j) {

			// alphaが小さい場合、simpsonによる積分が適さない
			double alpha = random->uniform(0.1, 1.0);
			glm::dvec3 Ng(0.0, 0.0, 1.0);
			glm::dvec3 wo = LambertianSampler::sample(random, Ng);

			MicrofacetCoupledDielectricsMaterial m;
			m.Ng = Ng;
			m.alpha = alpha;

			OnlineMean<double> mean;

			for (int i = 0; i < 500000; ++i) {
				glm::dvec3 wi = m.sample(random, wo);
				glm::dvec3 bxdf = m.bxdf(wo, wi);
				double pdf = m.pdf(wo, wi);
				double cosTheta = glm::dot(m.Ng, wi);

				REQUIRE(std::abs(bxdf.x - bxdf.y) < 1.0e-6);
				REQUIRE(std::abs(bxdf.y - bxdf.z) < 1.0e-6);

				glm::dvec3 value;
				if (glm::any(glm::greaterThanEqual(bxdf, glm::dvec3(1.0e-6f)))) {
					value = bxdf * cosTheta / pdf;
				}
				else {
					value = glm::dvec3(0.0);
				}

				mean.addSample(value.x);
			}
			double result = mean.mean();

			CAPTURE(alpha);
			CAPTURE(glm::dot(Ng, wo));
			REQUIRE(std::abs(result - 1.0) < 1.0e-2);
			printf("MC coupled dielectrics %.10f\n", result);
		}
	}
}

// 

TEST_CASE("ArbitraryBRDFSpace", "[ArbitraryBRDFSpace]") {
	using namespace rt;


	SECTION("ArbitraryBRDFSpace") {
		rt::Xor64 *random = new rt::Xor64();
		for (int j = 0; j < 1000000; ++j) {
			auto zAxis = sample_on_unit_sphere(random);
			ArbitraryBRDFSpace space(zAxis);

			REQUIRE(glm::abs(glm::dot(space.xaxis, space.yaxis)) < 1.0e-15);
			REQUIRE(glm::abs(glm::dot(space.yaxis, space.zaxis)) < 1.0e-15);
			REQUIRE(glm::abs(glm::dot(space.zaxis, space.xaxis)) < 1.0e-15);

			glm::dvec3 maybe_zaxis = glm::cross(space.xaxis, space.yaxis);
			for (int j = 0; j < 3; ++j) {
				REQUIRE(glm::abs(space.zaxis[j] - maybe_zaxis[j]) < 1.0e-15);
			}

			auto anyvector = sample_on_unit_sphere(random);
			auto samevector = space.localToGlobal(space.globalToLocal(anyvector));

			for (int j = 0; j < 3; ++j) {
				REQUIRE(glm::abs(anyvector[j] - samevector[j]) < 1.0e-15);
			}
		}
	}
}


TEST_CASE("Velvet", "[Velvet]") {
	using namespace rt;
	rt::XoroshiroPlus128 random;
	SECTION("D Normalization") {
		for (int j = 0; j < 128; ++j) {
			// alphaが小さい場合、simpsonによる積分が適さない
			double alpha = random.uniform(0.1, 1.0);

			double result = hemisphere_composite_simpson<double>([&](double theta, double phi) {
				glm::dvec3 wi = rt::polar_to_cartesian((double)theta, (double)phi);
				glm::dvec3 Ng(0.0, 0.0, 1.0);
				glm::dvec3 sample_m = rt::polar_to_cartesian((double)theta, (double)phi);
				double cosTheta = glm::dot(Ng, sample_m);
				double value = rt::velvet_D(Ng, sample_m, alpha) * cosTheta;
				return (double)value;
			}, 64);
			CAPTURE(alpha);
			CAPTURE(result);
			REQUIRE(glm::abs(result - 1.0) < 1.0e-5);
		}
	}
	SECTION("visible normal normalization") {
		glm::dvec3 Ng(0, 0, 1);
		for (int j = 0; j < 32; ++j) {
			double alpha = random.uniform(0.5, 1.0);
			glm::dvec3 wo = LambertianSampler::sample(&random, Ng);
			double cosThetaO = wo.z;

			double result = hemisphere_composite_simpson<double>([&](double theta, double phi) {
				glm::dvec3 h = rt::polar_to_cartesian((double)theta, (double)phi);
				double D = velvet_D(Ng, h, alpha);
				double G = velvet_G1(cosThetaO, alpha);
				double value = G * std::max(glm::dot(wo, h), 0.0) * D / cosThetaO;
				return value;
			}, 1000);

			// フィッティングなので、だいたいしか合わない
			CAPTURE(alpha);
			CAPTURE(result);
			REQUIRE(glm::abs(result - 1.0) < 0.05);
			// printf("%f, a = %f, cosTheta = %f\n", result, alpha, cosThetaO);
		}
	}
}

int main(int argc, char* const argv[])
{
#if 1
	// テストを指定する場合
	char* custom_argv[] = {
		"",
		//"[microfacet sampling]"
		//"[ArbitraryBRDFSpace]"
		//"[microfacet]"
		"[Velvet]"
	};
	Catch::Session().run(sizeof(custom_argv) / sizeof(custom_argv[0]), custom_argv);
#else
	// 全部やる場合
	char* custom_argv[] = {
		"",
	};
	Catch::Session session;
	session.run(sizeof(custom_argv) / sizeof(custom_argv[0]), custom_argv);
#endif

	std::cin.get();
	return 0;
}
