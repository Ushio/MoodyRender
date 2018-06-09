#pragma once

#include <glm/glm.hpp>
#include <glm/ext.hpp>
#include "adaptive_simpson.hpp"
#include "composite_simpson.hpp"

namespace rt {
	// https://ja.wikipedia.org/wiki/%E6%A5%B5%E5%BA%A7%E6%A8%99%E7%B3%BB
	// phi (0 ~ 2pi)
	// theta (0 ~ pi): 全球
	// theta (0 ~ pi/2): 半球

	// (theta, phi) -> value
	template <class Real>
	Real hemisphere_adaptive_simpson(std::function<Real(Real, Real)> f, Real eps) {
		std::function<Real(Real)> f_phi = [&](Real phi) {
			std::function<Real(Real)> f_theta = [&](Real theta) {
				double jacobian = std::sin(theta);
				return f(theta, phi) * jacobian;
			};
			return rt::adaptive_simpson(rt::SimpsonRange<Real>(f_theta, Real(0.0), glm::pi<Real>() * Real(0.5)), eps);
		};
		return rt::adaptive_simpson(rt::SimpsonRange<Real>(f_phi, Real(0.0), Real(2.0) * glm::pi<Real>()), eps);
	}

	// (theta, phi) -> value
	template <class Real>
	Real sphere_adaptive_simpson(std::function<Real(Real, Real)> f, Real eps) {
		std::function<Real(Real)> f_phi = [&](Real phi) {
			std::function<Real(Real)> f_theta = [&](Real theta) {
				double jacobian = std::sin(theta);
				return f(theta, phi) * jacobian;
			};
			return rt::adaptive_simpson(rt::SimpsonRange<Real>(f_theta, Real(0.0), glm::pi<Real>() * Real(1.0)), eps);
		};
		return rt::adaptive_simpson(rt::SimpsonRange<Real>(f_phi, Real(0.0), Real(2.0) * glm::pi<Real>()), eps);
	}

	template <class Real>
	Real hemisphere_composite_simpson(std::function<Real(Real, Real)> f, int n) {
		return rt::composite_simpson<Real>([&](double phi) {
			return rt::composite_simpson<Real>([&](double theta) {
				double jacobian = std::sin(theta);
				return f(theta, phi) * jacobian;
			}, n, Real(0.0), glm::pi<Real>() * Real(0.5));
		}, n, Real(0.0), Real(2.0) * glm::pi<Real>());
	}
}