#pragma once

#include <cmath>
#include <glm/glm.hpp>
#include <glm/ext.hpp>

namespace rt {
	inline glm::dvec3 polar_to_cartesian(double theta, double phi) {
		double sinTheta = std::sin(theta);
		glm::dvec3 v = {
			sinTheta * std::cos(phi),
			sinTheta * std::sin(phi),
			cos(theta)
		};
		return v;
	};

	// zが上の座標系に移動する行列
	inline glm::mat3 to_bxdf_basis_transform(const glm::dvec3 &n) {
		glm::dvec3 xaxis;
		glm::dvec3 zaxis = n;
		glm::dvec3 yaxis;
		if (0.999 < glm::abs(zaxis.z)) {
			xaxis = glm::normalize(glm::cross(glm::dvec3(0.0, 1.0, 0.0), zaxis));
		}
		else {
			xaxis = glm::normalize(glm::cross(glm::dvec3(0.0, 0.0, 1.0), zaxis));
		}
		yaxis = glm::cross(zaxis, xaxis);
		return glm::transpose(glm::mat3(xaxis, yaxis, zaxis));
	}
	inline glm::dvec3 from_bxdf(const glm::dvec3 &n, const glm::dvec3 &bxdf_dir) {
		return glm::transpose(to_bxdf_basis_transform(n)) * bxdf_dir;
	}

}