#pragma once

#include <cmath>
#include <glm/glm.hpp>
#include <glm/ext.hpp>

namespace rt {
	// この辺も整理したい
	inline glm::vec3 polar_to_cartesian(float theta, float phi) {
		float sinTheta = std::sin(theta);
		glm::vec3 v = {
			sinTheta * std::cos(phi),
			sinTheta * std::sin(phi),
			cos(theta)
		};
		return v;
	};

	// zが上の座標系に移動する行列
	inline glm::mat3 to_bxdf_basis_transform(const glm::vec3 &n) {
		glm::vec3 xaxis;
		glm::vec3 zaxis = n;
		glm::vec3 yaxis;
		if (0.999f < glm::abs(zaxis.z)) {
			xaxis = glm::normalize(glm::cross(glm::vec3(0.0f, 1.0f, 0.0f), zaxis));
		}
		else {
			xaxis = glm::normalize(glm::cross(glm::vec3(0.0f, 0.0f, 1.0f), zaxis));
		}
		yaxis = glm::cross(zaxis, xaxis);
		return glm::transpose(glm::mat3(xaxis, yaxis, zaxis));
	}
	inline glm::vec3 from_bxdf(const glm::vec3 &n, const glm::vec3 &bxdf_dir) {
		return glm::transpose(to_bxdf_basis_transform(n)) * bxdf_dir;
	}

}