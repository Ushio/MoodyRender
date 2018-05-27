#pragma once

#include <glm/glm.hpp>
#include <glm/ext.hpp>

#include "plane.hpp"
#include "peseudo_random.hpp"


namespace rt {
	inline glm::vec2 uniform_in_unit_circle(PeseudoRandom *random) {
		glm::vec2 d;
		float sq = 0.0f;
		do {
			d.x = random->uniform(-1.0f, 1.0f);
			d.y = random->uniform(-1.0f, 1.0f);
			sq = glm::length2(d);
		} while (1.0f < sq);
		return d;
	}

	struct CameraSetting {
		float fovy = glm::radians(45.0f);

		glm::vec3 eye = glm::vec3(0.0f, 0.0f, 1.0f);
		glm::vec3 lookat = glm::vec3(0.0f, 0.0f, 0.0f);
		glm::vec3 up = glm::vec3(0.0f, 1.0f, 0.0f);

		int imageWidth = 4;
		int imageHeight = 3;

		float lensRadius = 0.1f;
		float focasDistance = 3.0f;

		// TODO 整理
		float tanThetaH() const {
			return std::tan(fovy * 0.5f);
		}
		float distanceS() const {
			return imageHeight / (2.0f * tanThetaH());
		}
		float tanThetaV() const {
			return (float)imageWidth * 0.5f / distanceS();
		}

		float widthV() const {
			// 相似関係より
			return imageWidth * focasDistance / distanceS();
		}
		float heightV() const {
			return imageHeight * focasDistance / distanceS();
		}
	};

	class Camera {
	public:
		Camera() {}
		Camera(const CameraSetting &setting) :_setting(setting) {
			_view = glm::lookAt(_setting.eye, _setting.lookat, _setting.up);
			_viewInverse = glm::inverse(_view);

			glm::mat3 r(_viewInverse);
			_right = r[0];
			_up = r[1];
			_back = r[2];

			_planeV = plane_from(front(), origin() + front() * setting.focasDistance);
		}

		glm::mat4 view() const {
			return _view;
		}
		glm::mat4 viewInverse() const {
			return _viewInverse;
		}
		glm::vec3 origin() const {
			return _viewInverse[3];
		}
		glm::vec3 front() const {
			return -_back;
		}
		glm::vec3 back() const {
			return _back;
		}

		glm::vec3 right() const {
			return _right;
		}
		glm::vec3 left() const {
			return -_right;
		}

		glm::vec3 up() const {
			return _up;
		}
		glm::vec3 down() const {
			return -_up;
		}

		int imageWidth() const {
			return _setting.imageWidth;
		}
		int imageHeight() const {
			return _setting.imageHeight;
		}

		CameraSetting setting() const {
			return _setting;
		}

		// ちょっと下と処理が重複
		glm::vec3 sampleLens(PeseudoRandom *random) const {
			glm::vec2 sample = uniform_in_unit_circle(random);
			float r = setting().lensRadius;
			return origin() + r * right() * sample.x + r * down() * sample.y;
		}

		void sampleRay(PeseudoRandom *random, int x, int y, glm::vec3 *o, glm::vec3 *d) const {
			glm::vec2 sample = uniform_in_unit_circle(random);
			float r = setting().lensRadius;
			glm::vec3 sampleLens = origin() + r * right() * sample.x + r * down() * sample.y;

			auto focalPlaneCenter = origin() + front() * setting().focasDistance;

			auto width = setting().widthV();
			auto height = setting().heightV();

			glm::vec3 LT = focalPlaneCenter
				+ left() * width * 0.5f
				+ up() * height * 0.5f;

			float stepPixel = width / setting().imageWidth;

			glm::vec3 PixelLT = LT + stepPixel * right() * (float)x + stepPixel * down() * (float)y;
			glm::vec3 sampleFocalPlane = PixelLT
				+ right() * random->uniform() * stepPixel
				+ down() * random->uniform() * stepPixel;

			*o = sampleLens;
			*d = glm::normalize(sampleFocalPlane - sampleLens);
		}
		float lensPDF() const {
			float r = _setting.lensRadius;
			return 1.0f / (glm::pi<float>() * r * r);
		}

		Plane planeV() const {
			return _planeV;
		}

		bool findPixel(const glm::vec3 &o, const glm::vec3 &d, int *x, int *y) const {
			// 逆向きの光線は無視
			if (0.0f < glm::dot(d, front())) {
				return false;
			}

			float tmin = std::numeric_limits<float>::max();
			if (intersect_ray_plane(o, d, planeV(), &tmin)) {
				glm::vec3 Vp = o + tmin * d;

				auto focalPlaneCenter = origin() + front() * setting().focasDistance;
				auto width = setting().widthV();
				auto height = setting().heightV();

				glm::vec3 LT = focalPlaneCenter
					+ left() * width * 0.5f
					+ up() * height * 0.5f;

				glm::vec3 dir = Vp - LT;
				float Vx = glm::dot(dir, right());
				float Vy = glm::dot(dir, down());

				float stepPixel = width / setting().imageWidth;
				int at_x = (int)floor(Vx / stepPixel);
				int at_y = (int)floor(Vy / stepPixel);
				if (0 <= at_x && at_x < setting().imageWidth) {
					if (0 <= at_y && at_y < setting().imageHeight) {
						*x = at_x;
						*y = at_y;
						return true;
					}
				}
			}
			return false;
		}

		float Wi(const glm::vec3 &x0, const glm::vec3 &x1, const glm::vec3 &n1) const {
			glm::vec3 x1_to_x0 = x0 - x1;
			glm::vec3 x0_to_x1 = -x1_to_x0;
			float G = glm::dot(front(), glm::normalize(x0_to_x1)) * glm::dot(n1, glm::normalize(x1_to_x0)) / glm::length2(x0_to_x1);
			float c = glm::dot(n1, x1_to_x0) / std::pow(glm::dot(front(), x0_to_x1), 3.0f);
			float sigma = imageWidth() * imageHeight() / (4.0f * setting().tanThetaV() * setting().tanThetaH());
			return 1.0f / G * c * sigma * lensPDF() * 1.0f;
		}

		CameraSetting _setting;
		glm::mat4 _view;
		glm::mat4 _viewInverse;

		glm::vec3 _right;
		glm::vec3 _up;
		glm::vec3 _back;

		Plane _planeV;
	};
}
