#pragma once

#include <glm/glm.hpp>
#include <glm/ext.hpp>
#include "peseudo_random.hpp"

namespace rt {
	inline glm::dvec3 sample_on_unit_sphere(PeseudoRandom *random) {
		double x1;
		double x2;
		double S;
		do {
			x1 = random->uniform(-1.0, 1.0);
			x2 = random->uniform(-1.0, 1.0);
			S = x1 * x1 + x2 * x2;
		} while (S >= 1.0);

		double two_sqrt_one_minus_s = 2.0 * std::sqrt(std::max(1.0 - S, 0.0));
		return glm::dvec3(
			x1 * two_sqrt_one_minus_s,
			x2 * two_sqrt_one_minus_s,
			1.0 - 2.0 * S);
	}

	// z up
	inline glm::dvec3 sample_on_unit_hemisphere(PeseudoRandom *random) {
		double x1;
		double x2;
		double S;
		do {
			x1 = random->uniform(-1.0, 1.0);
			x2 = random->uniform(-1.0, 1.0);
			S = x1 * x1 + x2 * x2;
		} while (S >= 1.0);

		double c = std::sqrt(2.0 - S);
		return glm::dvec3(
			x1 * c,
			x2 * c,
			1.0 - S);
	}
}