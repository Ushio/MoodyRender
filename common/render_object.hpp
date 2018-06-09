#pragma once

#include <vector>

#include <glm/glm.hpp>
#include <glm/ext.hpp>
// #include <mapbox/variant.hpp>
#include <strict_variant/variant.hpp>

#include "camera.hpp"
#include "material.hpp"

namespace rt {
	class Geometry {
	public:
		struct Point {
			glm::vec3 P;
		};
		struct Primitive {
			glm::ivec3 indices;
			Material material;
		};
		std::vector<Point> points;
		std::vector<Primitive> primitives;
	};

	class Scene {
	public:
		std::vector<Geometry> geometries;
		rt::Camera camera;
	};
}
