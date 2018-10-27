#pragma once
#include "geometry.hpp"
#include "peseudo_random.hpp"
#include "spherical_triangle_sampler.hpp"

namespace rt {
	struct Rectangle {
		Rectangle() {}
		Rectangle(glm::dvec3 s, glm::dvec3 ex, glm::dvec3 ey)
			:_s(s)
			, _ex(ex)
			, _ey(ey) {
			_exLength = glm::length(ex);
			_eyLength = glm::length(ey);
			_x = ex / _exLength;
			_y = ey / _eyLength;
			_z = glm::cross(_x, _y);
		}
		glm::dvec3 sample(double u, double v) const {
			return _s + _ex * u + _ey * v;
		}
		double area() const {
			return _exLength * _eyLength;
		}
		glm::dvec3 normal() const {
			return _z;
		}
		glm::dvec3 s() const {
			return _s;
		}
		glm::dvec3 ex() const {
			return _ex;
		}
		glm::dvec3 ey() const {
			return _ey;
		}
		glm::dvec3 x() const {
			return _x;
		}
		glm::dvec3 y() const {
			return _y;
		}
		glm::dvec3 z() const {
			return _z;
		}
		double exLength() const {
			return _exLength;
		}
		double eyLength() const {
			return _eyLength;
		}
		glm::dvec3 _s;
		glm::dvec3 _ex;
		glm::dvec3 _ey;
		double _exLength = 0.0;
		double _eyLength = 0.0;
		glm::dvec3 _x;
		glm::dvec3 _y;
		glm::dvec3 _z;
	};

	class SphericalRectangleSamplerCoordinate_Optimized {
	public:
		SphericalRectangleSamplerCoordinate_Optimized(const Rectangle &rectangle, const glm::dvec3 &o) :_rectangle(rectangle) {
			_o = o;

			glm::dvec3 d = rectangle.s() - o;
			_x0 = glm::dot(d, rectangle.x());
			_y0 = glm::dot(d, rectangle.y());
			_z0 = glm::dot(d, rectangle.z());

			double exLen = rectangle.exLength();
			double eyLen = rectangle.eyLength();
			_x1 = _x0 + exLen;
			_y1 = _y0 + eyLen;

			// z flip
			if (_z0 > 0.0) {
				_z0 *= -1.0;
				_rectangle._z *= -1.0;
			}
			/*
			_v[0][0] = glm::vec3(_x0, _y0, _z0);
			_v[0][1] = glm::vec3(_x0, _y1, _z0);
			_v[1][0] = glm::vec3(_x1, _y0, _z0);
			_v[1][1] = glm::vec3(_x1, _y1, _z0);
			*/
			/*
			a.y * b.z - b.y * a.z,
			a.z * b.x - b.z * a.x,
			a.x * b.y - b.x * a.y
			*/
			// 外積を展開
			//double z0_exLen = _z0 * exLen;
			//double z0_exLen2 = z0_exLen * z0_exLen;
			//double z0_eyLen = _z0 * eyLen;
			//double z0_eyLen2 = z0_eyLen * z0_eyLen;
			//_n[0] = glm::dvec3(
			//	0.0,
			//	z0_exLen,
			//	-exLen * _y0
			//);
			//_n[1] = glm::dvec3(
			//	-z0_eyLen,
			//	0.0,
			//	_x1 * eyLen
			//);
			//_n[2] = glm::dvec3(
			//	0.0,
			//	-z0_exLen,
			//	_y1 * exLen
			//);
			//_n[3] = glm::dvec3(
			//	z0_eyLen,
			//	0.0,
			//	-_x0 * eyLen
			//);
			// 必要なのは、zだけであるので、正規化はzだけ
			//_n[0].z /= std::sqrt(z0_exLen2 + _n[0].z * _n[0].z);
			//_n[1].z /= std::sqrt(z0_eyLen2 + _n[1].z * _n[1].z);
			//_n[2].z /= std::sqrt(z0_exLen2 + _n[2].z * _n[2].z);
			//_n[3].z /= std::sqrt(z0_eyLen2 + _n[3].z * _n[3].z);
			// 定数倍なのだから、係数を掛ける必要が無い
			_z0z0 = _z0 * _z0;
			_y1y1 = _y1 * _y1;
			_n[0] = glm::dvec3(
				0.0,
				_z0,
				-_y0
			);
			_n[1] = glm::dvec3(
				-_z0,
				0.0,
				_x1
			);
			_n[2] = glm::dvec3(
				0.0,
				-_z0,
				_y1
			);
			_n[3] = glm::dvec3(
				_z0,
				0.0,
				-_x0
			);

			/*
			a.y * b.z - b.y * a.z,
			a.z * b.x - b.z * a.x,
			a.x * b.y - b.x * a.y

			a.y * b.z - 0.0 * a.z,
			a.z * b.x - b.z * 0.0,
			0.0 * 0.0 - b.x * a.y
			*/

			// 必要なのは、zだけであるので、正規化はzだけ
			_y0y0 = _y0 * _y0;
			_n[0].z /= std::sqrt(_z0z0 + _y0y0);
			_n[1].z /= std::sqrt(_z0z0 + _x1 * _x1);
			_n[2].z /= std::sqrt(_z0z0 + _y1y1);
			_n[3].z /= std::sqrt(_z0z0 + _x0 * _x0);

			/*
			_n[0].x == 0
			_n[1].y == 0
			_n[2].x == 0
			_n[3].y == 0
			であるため、z成分だけで良い
			*/
			_g[0] = std::acos(-_n[0].z * _n[1].z);
			_g[1] = std::acos(-_n[1].z * _n[2].z);
			_g[2] = std::acos(-_n[2].z * _n[3].z);
			_g[3] = std::acos(-_n[3].z * _n[0].z);

			_sr = _g[0] + _g[1] + _g[2] + _g[3] - glm::two_pi<double>();
		}

		double solidAngle() const {
			return _sr;
		}

		glm::dvec3 sample(double u, double v) const {
			double AQ = _sr;
			double phi_u = u * AQ - _g[2] - _g[3] + glm::two_pi<double>();

			auto safeSqrt = [](double x) {
				return std::sqrt(std::max(x, 0.0));
			};

			double b0 = _n[0].z;
			double b1 = _n[2].z;
			double fu = (std::cos(phi_u) * b0 - b1) / std::sin(phi_u);
			double cu = std::copysign(1.0, fu) / std::sqrt(fu * fu + b0 * b0);
			double xu = -cu * _z0 / safeSqrt(1.0 - cu * cu);

			double d = std::sqrt(xu * xu + _z0z0);
			double d2 = d * d;
			double h0 = _y0 / std::sqrt(d2 + _y0y0);
			double h1 = _y1 / std::sqrt(d2 + _y1y1);
			double hv = glm::mix(h0, h1, v);
			double yv = hv * d / safeSqrt(1.0 - hv * hv);
			return _o + xu * _rectangle.x() + yv * _rectangle.y() + _z0 * _rectangle.z();
		}

		Rectangle _rectangle;

		glm::dvec3 _o;

		double _x0;
		double _x1;
		double _y0;
		double _y1;
		double _z0;

		double _z0z0;
		double _y0y0;
		double _y1y1;

		double _sr;

		glm::dvec3 _n[4];
		double _g[4];
	};
	using SphericalRectangleSamplerCoordinate = SphericalRectangleSamplerCoordinate_Optimized;

	class IDirectSampler {
	public:
		// サンプリングが可能かどうか can_sample == falseなら pdf = 0である
		virtual bool can_sample(glm::dvec3 o) const = 0;
		virtual double pdf_area(glm::dvec3 o, glm::dvec3 p) const = 0;

		// can_sample == false のときは呼ばないことにする
		virtual void sample(PeseudoRandom *random, glm::dvec3 o, glm::dvec3 *p, glm::dvec3 *n, glm::dvec3 *Le, double *pdf_area) const = 0;

		// for light selection
		virtual glm::dvec3 center() const = 0;
		virtual double Lavg_mul_area() const = 0;
	};

	/*
	表面から
	*/
	class LightSelectorHeuristic {
	public:
		LightSelectorHeuristic(glm::dvec3 o, std::vector<IDirectSampler *>::const_iterator beg, std::vector<IDirectSampler *>::const_iterator end, std::vector<double> &importance_storage)
			: _beg(beg), _end(end), _importances(importance_storage) {
			
			_importances.clear();
			for (auto it = _beg; it != _end; ++it) {
				double importance;
				if ((*it)->can_sample(o)) {
					_can_sample = true;
					importance = (*it)->Lavg_mul_area() / glm::distance2((*it)->center(), o);
					_importance_sum += importance;
				}
				else {
					importance = 0.0;
				}
				_importances.push_back(importance);
			}
		}
		bool can_sample() const {
			return _can_sample;
		}
		const IDirectSampler *choice(PeseudoRandom *random, double *p_choice) const {
			// どうせこの手法はスケールしない。
			// なので割り切ってリニアサーチにする
			double s = random->uniform(0.0, _importance_sum);
			for (int i = 0; i < _importances.size(); ++i) {
				if (s < _importances[i]) {
					*p_choice = p(i);
					return _beg[i];
				}
				s -= _importances[i];
			}
			int index = std::distance(_beg, _end) - 1;

			*p_choice = p(index);
			return _beg[index];
		}
		double p(int i) const {
			return _importances[i] / _importance_sum;
		}
		double p(const IDirectSampler *sampler) const {
			for (auto it = _beg; it != _end; ++it) {
				if ((*it) == sampler) {
					return p(std::distance(_beg, it));
				}
			}
			return 0.0;
		}
	private:
		std::vector<IDirectSampler *>::const_iterator _beg;
		std::vector<IDirectSampler *>::const_iterator _end;
		bool _can_sample = false;
		std::vector<double> &_importances;
		double _importance_sum = 0.0;
	};
	class LightSelectorUniform {
	public:
		LightSelectorUniform(glm::dvec3 o, std::vector<IDirectSampler *>::const_iterator beg, std::vector<IDirectSampler *>::const_iterator end, std::vector<double> &importance_storage)
			: _beg(beg), _end(end) {
		}
		bool can_sample() const {
			return true;
		}
		const IDirectSampler *choice(PeseudoRandom *random, double *p_choice) const {
			int n = std::distance(_beg, _end);
			double s = random->uniform(0.0, n);
			int index = std::floor(s);
			index = std::min(index, n - 1);
			*p_choice = p(index);
			return _beg[index];
		}
		double p(int i) const {
			int n = std::distance(_beg, _end);
			return 1.0 / n;
		}
		double p(const IDirectSampler *sampler) const {
			for (auto it = _beg; it != _end; ++it) {
				if ((*it) == sampler) {
					return p(std::distance(_beg, it));
				}
			}
			return 0.0;
		}
	private:
		std::vector<IDirectSampler *>::const_iterator _beg;
		std::vector<IDirectSampler *>::const_iterator _end;
	};
	using LightSelector = LightSelectorHeuristic;
	// using LightSelector = LightSelectorUniform;

	class TriangleAreaSampler : public IDirectSampler {
	public:
		TriangleAreaSampler(glm::dvec3 a, glm::dvec3 b, glm::dvec3 c, bool doubleSided, glm::dvec3 Le)
			:_a(a)
			, _b(b)
			, _c(c)
			, _doubleSided(doubleSided)
			, _Le(Le) {
			_Lavg = (_Le.x + _Le.y + _Le.z) / 3.0;
			_n = triangleNormal(_a, _b, _c);
			_center = (_a + _b + _c) / 3.0;
			_area = triangleArea(_a, _b, _c);
			_one_over_area = 1.0 / _area;

			_Lavg_mul_area = _Lavg * _area;
		}

		virtual double pdf_area(glm::dvec3 o, glm::dvec3 p) const override {
			if (can_sample(o) == false) {
				return 0.0;
			}
			return _one_over_area;
		}
		virtual bool can_sample(glm::dvec3 o) const override {
			// サンプル面からの符号付き距離
			// sd > 0.0 なら表
			double sd = glm::dot(o - _a, _n);
			const double kEps = 1.0e-6;

			if (_doubleSided) {
				// 面に対して水平でない : |sd| > kEps
				return std::abs(sd) > kEps;
			}
			else {
				// 面に対して水平でない : |sd| > kEps
				// 表                   : sd > 0
				// したがって、
				// sd > kEps
				return sd > kEps;
			}
		}

		virtual void sample(PeseudoRandom *random, glm::dvec3 o, glm::dvec3 *p, glm::dvec3 *n, glm::dvec3 *Le, double *pdf_area) const override {
			*p = uniform_on_triangle(random->uniform(), random->uniform()).evaluate(_a, _b, _c);

			glm::dvec3 d = *p - o;
			bool backfacing = 0.0 < glm::dot(d, _n);

			if (_doubleSided) {
				*n = backfacing ? -_n : _n;
				*Le = _Le;
			}
			else {
				*n = _n;
				*Le = backfacing ? glm::dvec3(0.0) : _Le;
			}
			*pdf_area = _one_over_area;
		}

		virtual glm::dvec3 center() const override { return _center; }
		virtual double Lavg_mul_area() const override { return _Lavg_mul_area; }
	private:
		glm::dvec3 _Le;
		double _Lavg = 0;
		glm::dvec3 _a, _b, _c;
		glm::dvec3 _n;
		glm::dvec3 _center;
		bool _doubleSided = false;
		double _area = 0.0;
		double _one_over_area = 0.0;

		double _Lavg_mul_area = 0.0;
	};
	class SphericalRectangleSampler : public IDirectSampler {
	public:
		SphericalRectangleSampler(glm::dvec3 s, glm::dvec3 ex, glm::dvec3 ey, bool doubleSided, glm::dvec3 Le)
			:_doubleSided(doubleSided)
			, _Le(Le)
			, _q(s, ex, ey)
		{
			_Lavg = (_Le.x + _Le.y + _Le.z) / 3.0;
			_center = s + ex * 0.5 + ey * 0.5;

			_Lavg_mul_area = _Lavg * _q.area();
		}
		virtual bool can_sample(glm::dvec3 o) const override {
			// サンプル面からの符号付き距離
			// sd > 0.0 なら表
			double sd = glm::dot(o - _q.s(), _q.normal());
			const double kEps = 1.0e-6;

			if (_doubleSided) {
				// 面に対して水平でない : |sd| > kEps
				return std::abs(sd) > kEps;
			}
			else {
				// 面に対して水平でない : |sd| > kEps
				// 表                   : sd > 0
				// したがって、
				// sd > kEps
				return sd > kEps;
			}
		}
		virtual void sample(PeseudoRandom *random, glm::dvec3 o, glm::dvec3 *p, glm::dvec3 *n, glm::dvec3 *Le, double *pdf_area) const
		{
			SphericalRectangleSamplerCoordinate sampler(_q, o);
			*p = sampler.sample(random->uniform(), random->uniform());
			glm::dvec3 d = *p - o;
			bool backfacing = 0.0 < glm::dot(d, _q.normal());

			if (_doubleSided) {
				*n = backfacing ? -_q.normal() : _q.normal();
				*Le = _Le;
			}
			else {
				*n = _q.normal();
				*Le = backfacing ? glm::dvec3(0.0) : _Le;
			}

			double dLength2 = glm::length2(d);
			double cosTheta = glm::dot(-*n, d / std::sqrt(dLength2));
			double pw = 1.0 / sampler.solidAngle();
			*pdf_area = pw * cosTheta / dLength2;
		}
		virtual double pdf_area(glm::dvec3 o, glm::dvec3 p) const {
			if (can_sample(o) == false) {
				return 0.0;
			}

			glm::dvec3 d = p - o;
			bool backfacing = 0.0 < glm::dot(d, _q.normal());
			glm::dvec3 n = backfacing ? -_q.normal() : _q.normal();

			SphericalRectangleSamplerCoordinate sampler(_q, o);
			double dLength2 = glm::length2(d);
			double cosTheta = glm::dot(-n, d / std::sqrt(dLength2));
			double pw = 1.0 / sampler.solidAngle();
			return pw * cosTheta / glm::distance2(o, p);
		}

		virtual glm::dvec3 center() const override { return _center; }
		virtual double Lavg_mul_area() const override { return _Lavg_mul_area; }
	private:
		glm::dvec3 _Le;
		double _Lavg = 0;
		bool _doubleSided = false;
		Rectangle _q;
		glm::dvec3 _center;

		double _Lavg_mul_area = 0.0;
	};


	class SphericalTriangleDirectSampler : public IDirectSampler {
	public:
		SphericalTriangleDirectSampler(glm::dvec3 a, glm::dvec3 b, glm::dvec3 c, bool doubleSided, glm::dvec3 Le)
			:_a(a)
			, _b(b)
			, _c(c)
			, _doubleSided(doubleSided)
			, _Le(Le) {
			_Lavg = (_Le.x + _Le.y + _Le.z) / 3.0;
			_n = triangleNormal(_a, _b, _c);
			_center = (_a + _b + _c) / 3.0;
			double area = triangleArea(_a, _b, _c);
			_Lavg_mul_area = _Lavg * area;
		}

		virtual double pdf_area(glm::dvec3 o, glm::dvec3 p) const override {
			if (can_sample(o) == false) {
				return 0.0;
			}
			
			glm::dvec3 d = p - o;
			bool backfacing = 0.0 < glm::dot(d, _n);
			glm::dvec3 n = backfacing ? -_n : _n;

			SphericalTriangleSampler sampler(_a, _b, _c, _n, o);
			double dLength2 = glm::length2(d);
			double cosTheta = glm::dot(-n, d / std::sqrt(dLength2));
			double pw = 1.0 / sampler.solidAngle();
			return pw * cosTheta / dLength2;
		}
		virtual bool can_sample(glm::dvec3 o) const override {
			// サンプル面からの符号付き距離
			// sd > 0.0 なら表
			double sd = glm::dot(o - _a, _n);
			const double kEps = 1.0e-6;

			if (_doubleSided) {
				// 面に対して水平でない : |sd| > kEps
				return std::abs(sd) > kEps;
			}
			else {
				// 面に対して水平でない : |sd| > kEps
				// 表                   : sd > 0
				// したがって、
				// sd > kEps
				return sd > kEps;
			}
		}

		virtual void sample(PeseudoRandom *random, glm::dvec3 o, glm::dvec3 *p, glm::dvec3 *n, glm::dvec3 *Le, double *pdf_area) const override {
			SphericalTriangleSampler sampler(_a, _b, _c, _n, o);

			// dは単位ベクトル
			glm::dvec3 d;
			*p = sampler.sample(random->uniform(), random->uniform(), &d);

			bool backfacing = 0.0 < glm::dot(d, _n);

			if (_doubleSided) {
				*n = backfacing ? -_n : _n;
				*Le = _Le;
			}
			else {
				*n = _n;
				*Le = backfacing ? glm::dvec3(0.0) : _Le;
			}

			double dLength2 = glm::distance2(*p, o);
			double cosTheta = glm::dot(-*n, d);
			double pw = 1.0 / sampler.solidAngle();
			*pdf_area = pw * cosTheta / dLength2;
		}

		virtual glm::dvec3 center() const override { return _center; }
		virtual double Lavg_mul_area() const override { return _Lavg_mul_area; }
	private:
		glm::dvec3 _Le;
		double _Lavg = 0;
		glm::dvec3 _a, _b, _c;
		glm::dvec3 _n;
		glm::dvec3 _center;
		bool _doubleSided = false;
		double _Lavg_mul_area = 0.0;
	};
}