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
			float alpha = random->uniformf(0.1f, 1.0f);
			glm::vec3 Ng(0.0f, 0.0f, 1.0f);
			glm::vec3 wo = LambertianSampler::sample(random, Ng);

			LambertianMaterial m;
			m.Ng = Ng;
			m.Le = glm::vec3(0.0);
			m.R = glm::vec3(1.0);

			double result = hemisphere_composite_simpson<double>([&](double theta, double phi) {
				glm::vec3 wi = rt::polar_to_cartesian((float)theta, (float)phi);

				glm::vec3 brdf = bxdf_evaluate(m, wo, wi);
				float cosTheta = glm::dot(bxdf_Ng(m), wi);

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
		ofToDataPath("baked/albedo_specular_conductor.xml").c_str(),
		ofToDataPath("baked/albedo_specular_conductor_avg.xml").c_str());
	rt::CoupledBRDFDielectrics::load(
		ofToDataPath("baked/albedo_specular_dielectrics.xml").c_str(),
		ofToDataPath("baked/albedo_specular_dielectrics_avg.xml").c_str());
	//rt::CoupledBRDFConductor::load(ofToDataPath("baked/albedo_specular_conductor.xml").c_str(), ofToDataPath("baked/albedo_specular_conductor_avg.xml").c_str());
	//rt::CoupledBRDFDielectrics::load(ofToDataPath("baked/albedo_specular_dielectrics.xml").c_str(), ofToDataPath("baked/albedo_specular_dielectrics_avg.xml").c_str());

	//SECTION("hemisphere_composite_simpson") {
	//	double result = rt::hemisphere_composite_simpson<double>([](double theta, double phi) {
	//		return 1.0 / glm::pi<double>() * std::cos(theta);
	//	}, 100);
	//	REQUIRE(std::abs(result - 1.0) < 1.0e-8);
	//}

	//SECTION("beckmann normalization") {
	//	using namespace rt;

	//	rt::Xor64 random;
	//	for (int j = 0; j < 100; ++j) {
	//		//int sample = 0;

	//		// alphaが小さい場合、simpsonによる積分が適さない
	//		float alpha = random.uniformf(0.1f, 1.0f);

	//		double result = hemisphere_composite_simpson<double>([&](double theta, double phi) {
	//			glm::vec3 wi = rt::polar_to_cartesian((float)theta, (float)phi);
	//			glm::vec3 Ng(0.0f, 0.0f, 1.0f);
	//			glm::vec3 sample_m = rt::polar_to_cartesian((float)theta, (float)phi);
	//			float cosTheta = glm::dot(Ng, sample_m);
	//			float value = rt::D_Beckmann(Ng, sample_m, alpha) * cosTheta;
	//			return (double)value;
	//		}, 500);
	//		CAPTURE(result);
	//		CAPTURE(alpha);
	//		REQUIRE(std::abs(result - 1.0) < 1.0e-4);

	//	}
	//}

	//SECTION("Importance Sampling MicrofacetCoupledConductorMaterial") {
	//	using namespace rt;
	//	rt::Xor64 *random = new rt::Xor64();
	//	for (int j = 0; j < 100; ++j) {
	//		OnlineMean<double> mean;
	//		float alpha = random->uniformf(0.1f, 1.0f);
	//		glm::vec3 Ng(0.0f, 0.0f, 1.0f);

	//		for (int i = 0; i < 1000000; ++i) {
	//			//float theta = CoupledBRDFConductor::sampler().sampleTheta(alpha, random);
	//			//glm::vec3 sample = polar_to_cartesian(theta, random->uniformf(0.0f, glm::two_pi<float>()));
	//			//ArbitraryBRDFSpace space(Ng);
	//			//glm::vec3 wi = space.localToGlobal(sample);

	//			// glm::vec3 wi = LambertianSampler::sample(random, Ng);

	//			const CoupledBRDFSampler &sampler = CoupledBRDFConductor::sampler();
	//			// float pdf = (1.0f / glm::two_pi<float>()) * (2.0 * sampler.thetaSize(alpha) / glm::pi<float>()) * sampler.probability(alpha, theta) / std::sin(theta);
	//			// float pdf = LambertianSampler::pdf(wi, Ng);
	//			// float value = 1.0 / glm::pi<float>() * glm::dot(wi, Ng) / pdf;
	//			//float theta = random->uniformf(0.0f, glm::pi<float>() * 0.5f);
	//			//// float p = 1.0 / sampler.thetaSize(alpha);
	//			//float p = sampler.probability(alpha, theta);
	//			//float pdf = p * sampler.thetaSize(alpha) * (2.0 / glm::pi<float>());
	//			//float value = pdf / (2.0 / glm::pi<float>());


	//			glm::vec3 wi = LambertianSampler::sample(random, Ng);
	//			float theta = acos(glm::dot(wi, Ng));

	//			float value = (1.0f / glm::two_pi<float>()) * (sampler.thetaSize(alpha) * (2.0 / glm::pi<float>())) * sampler.probability(alpha, theta) / std::sin(theta);

	//			float sample = value / LambertianSampler::pdf(wi, Ng);
	//			if (std::isfinite(sample) == false) {
	//				printf("\n");
	//			}
	//			mean.addSample(sample);
	//		}

	//		float result = mean.mean();
	//		printf("%f\n", result);
	//	}
	//}
	SECTION("white furnance test MicrofacetCoupledConductorMaterial") {
		using namespace rt;

		rt::Xor64 *random = new rt::Xor64();
		for (int j = 0; j < 100; ++j) {

			// alphaが小さい場合、simpsonによる積分が適さない
			float alpha = random->uniformf(0.1f, 1.0f);
			glm::vec3 Ng(0.0f, 0.0f, 1.0f);
			glm::vec3 wo = LambertianSampler::sample(random, Ng);

			MicrofacetCoupledConductorMaterial m;
			m.Ng = Ng;
			m.alpha = alpha;
			m.useFresnel = false;

			double result = hemisphere_composite_simpson<double>([&](double theta, double phi) {
				glm::vec3 wi = rt::polar_to_cartesian((float)theta, (float)phi);

				glm::vec3 brdf = bxdf_evaluate(m, wo, wi);
				float cosTheta = glm::dot(bxdf_Ng(m), wi);

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

	SECTION("white furnance test MC MicrofacetCoupledConductorMaterial") {
		using namespace rt;

		rt::Xor64 *random = new rt::Xor64();
		for (int j = 0; j < 32; ++j) {

			// alphaが小さい場合、simpsonによる積分が適さない
			float alpha = random->uniformf(0.1f, 1.0f);
			glm::vec3 Ng(0.0f, 0.0f, 1.0f);
			glm::vec3 wo = LambertianSampler::sample(random, Ng);

			MicrofacetCoupledConductorMaterial m;
			m.Ng = Ng;
			m.alpha = alpha;
			m.useFresnel = false;

			OnlineMean<double> mean;

			for (int i = 0; i < 500000; ++i) {
				glm::vec3 wi = bxdf_sample(m, random, wo);
				glm::vec3 bxdf = bxdf_evaluate(m, wo, wi);
				float pdf = bxdf_pdf(m, wo, wi);
				float cosTheta = glm::dot(bxdf_Ng(m), wi);

				REQUIRE(std::abs(bxdf.x - bxdf.y) < 1.0e-6);
				REQUIRE(std::abs(bxdf.y - bxdf.z) < 1.0e-6);

				glm::vec3 value;
				if (glm::any(glm::greaterThanEqual(bxdf, glm::vec3(1.0e-6f)))) {
					value = bxdf * cosTheta / pdf;
				}
				else {
					value = glm::vec3(0.0f);
				}

				mean.addSample(value.x);
			}
			double result = mean.mean();

			CAPTURE(alpha);
			CAPTURE(glm::dot(Ng, wo));
			REQUIRE(std::abs(result - 1.0) < 1.0e-2);
			// printf("%f\n", result);
		}
	}


}



int main(int argc, char* const argv[])
{
#if 1
	// テストを指定する場合
	char* custom_argv[] = {
		"",
		"[microfacet]"
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
