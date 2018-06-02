#include "ofMain.h"

#define CATCH_CONFIG_RUNNER
#include <catch.hpp>

#include <set>

#include "online.hpp"
#include "peseudo_random.hpp"
#include "composite_simpson.hpp"
#include "adaptive_simpson.hpp"

TEST_CASE("online", "[online]") {
	SECTION("online") {
		rt::Xor random;
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

	SECTION("bad sample adaptive_simpson") {
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

		rt::Xor random;

		SECTION("composite_simpson") {
			for (int i = 0; i < 100000; ++i) {
				double a = random.uniform(0, 5.0);
				double b = random.uniform(0, 5.0);
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
			for (int i = 0; i < 100000; ++i) {
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


int main(int argc, char* const argv[])
{
#if 0
	// テストを指定する場合
	char* custom_argv[] = {
		"",
		"[factorial]"
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
