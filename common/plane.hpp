#pragma once

#include <glm/glm.hpp>
#include <glm/ext.hpp>

namespace rt {
	// Plane equation
	// d = dot(n, p) for a given point p on the plane
	struct Plane {
		glm::vec3 n;
		float d = 0.0f;
	};
	inline Plane plane_from(const glm::vec3 &n, const glm::vec3 &p) {
		Plane plane;
		plane.n = n;
		plane.d = glm::dot(plane.n, p);
		return plane;
	}

	inline bool intersect_ray_plane(const glm::vec3 &o, const glm::vec3 &d, const Plane &plane, float *tmin) {
		float eps = 1.0e-5f;
		auto denom = glm::dot(plane.n, d);
		if (std::fabs(denom) < eps) {
			return false;
		}
		auto this_tmin = (plane.d - glm::dot(plane.n, o)) / denom;
		if (*tmin < this_tmin) {
			return false;
		}
		*tmin = this_tmin;
		return true;
	}

}