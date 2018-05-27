#pragma once

#include <vector>

#include <glm/glm.hpp>
#include <glm/ext.hpp>
// #include <mapbox/variant.hpp>
#include <strict_variant/variant.hpp>

#include "camera.hpp"

namespace rt {
	class LambertianMaterial {
	public:
		LambertianMaterial():Le(0.0f), R(1.0f) {}
		LambertianMaterial(glm::vec3 e, glm::vec3 r) : Le(e), R(r) {}
		glm::vec3 Le;
		glm::vec3 R;
		glm::vec3 Ng;

		bool isEmissive() const {
			return glm::any(glm::greaterThanEqual(Le, glm::vec3(glm::epsilon<float>())));
		}
	};
	class SpecularMaterial {
	public:
		glm::vec3 Ng;
	};

	typedef strict_variant::variant<LambertianMaterial, SpecularMaterial> Material;

	namespace MaterialVisitor {
		struct SetNg {
			SetNg(const glm::vec3 &Ng) : _Ng(Ng) {}
			template <class T>
			void operator()(T &m) {
				m.Ng = _Ng;
			}
			glm::vec3 _Ng;
		};
		struct GetNg {
			template <class T>
			glm::vec3 operator()(T &m) {
				return m.Ng;
			}
		};
	}

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
