#pragma once

#include <glm/glm.hpp>
#include <glm/ext.hpp>

#include "plane.hpp"
#include "peseudo_random.hpp"


namespace rt {
	inline glm::dvec2 uniform_in_unit_circle(PeseudoRandom *random) {
		glm::dvec2 d;
		double sq = 0.0;
		do {
			d.x = random->uniform(-1.0, 1.0);
			d.y = random->uniform(-1.0, 1.0);
			sq = glm::length2(d);
		} while (1.0 < sq);
		return d;
	}

	struct CameraSetting {
		double fovy = glm::radians(45.0);

		glm::dvec3 eye = glm::dvec3(0.0, 0.0, 1.0);
		glm::dvec3 lookat = glm::dvec3(0.0, 0.0, 0.0);
		glm::dvec3 up = glm::dvec3(0.0, 1.0, 0.0);

		int imageWidth = 4;
		int imageHeight = 3;

		double lensRadius = 0.1;
		double focasDistance = 3.0;

		// TODO 整理
		double tanThetaH() const {
			return std::tan(fovy * 0.5);
		}
		double distanceS() const {
			return imageHeight / (2.0 * tanThetaH());
		}
		double tanThetaV() const {
			return (double)imageWidth * 0.5 / distanceS();
		}

		double widthV() const {
			// 相似関係より
			return imageWidth * focasDistance / distanceS();
		}
		double heightV() const {
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

		glm::dmat4 view() const {
			return _view;
		}
		glm::dmat4 viewInverse() const {
			return _viewInverse;
		}
		glm::dvec3 origin() const {
			return _viewInverse[3];
		}
		glm::dvec3 front() const {
			return -_back;
		}
		glm::dvec3 back() const {
			return _back;
		}

		glm::dvec3 right() const {
			return _right;
		}
		glm::dvec3 left() const {
			return -_right;
		}

		glm::dvec3 up() const {
			return _up;
		}
		glm::dvec3 down() const {
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
		glm::dvec3 sampleLens(PeseudoRandom *random) const {
			glm::dvec2 sample = uniform_in_unit_circle(random);
			double r = setting().lensRadius;
			return origin() + r * right() * sample.x + r * down() * sample.y;
		}

		void sampleRay(PeseudoRandom *random, int x, int y, glm::dvec3 *o, glm::dvec3 *d) const {
			glm::dvec2 sample = uniform_in_unit_circle(random);
			double r = setting().lensRadius;
			glm::dvec3 sampleLens = origin() + r * right() * sample.x + r * down() * sample.y;

			auto focalPlaneCenter = origin() + front() * setting().focasDistance;

			auto width = setting().widthV();
			auto height = setting().heightV();

			glm::dvec3 LT = focalPlaneCenter
				+ left() * width * 0.5
				+ up() * height * 0.5;

			double stepPixel = width / setting().imageWidth;

			glm::dvec3 PixelLT = LT + stepPixel * right() * (double)x + stepPixel * down() * (double)y;
			glm::dvec3 sampleFocalPlane = PixelLT
				+ right() * random->uniform() * stepPixel
				+ down() * random->uniform() * stepPixel;

			*o = sampleLens;
			*d = glm::normalize(sampleFocalPlane - sampleLens);
		}
		double lensPDF() const {
			double r = _setting.lensRadius;
			return 1.0 / (glm::pi<double>() * r * r);
		}

		Plane planeV() const {
			return _planeV;
		}

		bool findPixel(const glm::dvec3 &o, const glm::dvec3 &d, int *x, int *y) const {
			// 逆向きの光線は無視
			if (0.0 < glm::dot(d, front())) {
				return false;
			}

			double tmin = std::numeric_limits<double>::max();
			if (intersect_ray_plane(o, d, planeV(), &tmin)) {
				glm::dvec3 Vp = o + tmin * d;

				auto focalPlaneCenter = origin() + front() * setting().focasDistance;
				auto width = setting().widthV();
				auto height = setting().heightV();

				glm::dvec3 LT = focalPlaneCenter
					+ left() * width * 0.5
					+ up() * height * 0.5;

				glm::dvec3 dir = Vp - LT;
				double Vx = glm::dot(dir, right());
				double Vy = glm::dot(dir, down());

				double stepPixel = width / setting().imageWidth;
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

		double Wi(const glm::dvec3 &x0, const glm::dvec3 &x1, const glm::dvec3 &n1) const {
			glm::dvec3 x1_to_x0 = x0 - x1;
			glm::dvec3 x0_to_x1 = -x1_to_x0;
			double G = glm::dot(front(), glm::normalize(x0_to_x1)) * glm::dot(n1, glm::normalize(x1_to_x0)) / glm::length2(x0_to_x1);
			double c = glm::dot(n1, x1_to_x0) / std::pow(glm::dot(front(), x0_to_x1), 3.0);
			double sigma = imageWidth() * imageHeight() / (4.0 * setting().tanThetaV() * setting().tanThetaH());
			return 1.0 / G * c * sigma * lensPDF() * 1.0;
		}

		CameraSetting _setting;
		glm::dmat4 _view;
		glm::dmat4 _viewInverse;

		glm::dvec3 _right;
		glm::dvec3 _up;
		glm::dvec3 _back;

		Plane _planeV;
	};
}
