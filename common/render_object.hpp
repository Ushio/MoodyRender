#pragma once

#include <vector>

#include <glm/glm.hpp>
#include <glm/ext.hpp>
#include <mapbox/variant.hpp>

#include "camera.hpp"

namespace rt {
	class LambertianMaterial {
	public:
		LambertianMaterial() {}
		LambertianMaterial(glm::vec3 e, glm::vec3 r) : Le(e), R(r) {}
		glm::vec3 Le;
		glm::vec3 R;
		glm::vec3 Ng;

		bool isEmissive() const {
			return glm::any(glm::greaterThanEqual(Le, glm::vec3(glm::epsilon<float>())));
		}
	};
	static const char *LambertianMaterialString = "LambertianMaterial";

	class SpecularMaterial {
	public:
		glm::vec3 Ng;
	};
	static const char *SpecularMaterialString = "SpecularMaterial";

	typedef mapbox::util::variant<LambertianMaterial, SpecularMaterial> Material;

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
