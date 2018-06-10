#pragma once

#include <glm/glm.hpp>
#include <glm/ext.hpp>

namespace rt {
	inline glm::vec3 triNg(const glm::vec3 &v0, const glm::vec3 &v1, const glm::vec3 &v2, bool isback) {
		glm::vec3 e1 = v1 - v0;
		glm::vec3 e2 = v2 - v0;
		return glm::normalize(isback ? glm::cross(e2, e1) : glm::cross(e1, e2));
	}
	inline double triArea(const glm::vec3 &p0, const glm::vec3 &p1, const glm::vec3 &p2) {
		auto va = p0 - p1;
		auto vb = p2 - p1;
		return glm::length(glm::cross(va, vb)) * 0.5;
	}
}