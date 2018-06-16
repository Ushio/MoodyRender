#pragma once

#include <glm/glm.hpp>
#include <glm/ext.hpp>

#include "peseudo_random.hpp"

namespace rt {
	inline glm::dvec3 triNg(const glm::dvec3 &v0, const glm::dvec3 &v1, const glm::dvec3 &v2, bool isback) {
		glm::dvec3 e1 = v1 - v0;
		glm::dvec3 e2 = v2 - v0;
		return glm::normalize(isback ? glm::cross(e2, e1) : glm::cross(e1, e2));
	}
	inline double triArea(const glm::dvec3 &p0, const glm::dvec3 &p1, const glm::dvec3 &p2) {
		auto va = p0 - p1;
		auto vb = p2 - p1;
		return glm::length(glm::cross(va, vb)) * 0.5;
	}

	template <class Real>
	struct TriangleSample {
		Real alpha = Real(0.0);
		Real beta = Real(0.0);

		template <class T>
		T evaluate(const T &A, const T &B, const T &C) const {
			return A * alpha + B * (Real(1.0) - beta) + C * (beta - alpha);
		}
	};

	// eps1, eps2: uniform 0~1
	template <class Real>
	inline TriangleSample<Real> uniform_on_triangle(Real eps1, Real eps2) {
		TriangleSample<Real> s;
		s.alpha = glm::min(eps1, eps2);
		s.beta = glm::max(eps1, eps2);
		return s;
	}
}