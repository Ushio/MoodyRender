#pragma once
#include <numeric>

namespace rt {
	// a * x + b
	template <class Real>
	class LinearTransform {
	public:
		LinearTransform() : _a(Real(1.0)), _b(Real(0.0)) {}

		LinearTransform(Real inputMin, Real inputMax, Real outputMin, Real outputMax) {
			_a = (outputMax - outputMin) / (inputMax - inputMin);
			_b = outputMin - _a * inputMin;
		}
		Real operator()(Real x) const {
			return std::fma(_a, x, _b);
		}
		bool isValid() const {
			return std::isfinite(_a) && std::isfinite(_b);
		}
		LinearTransform<double> inverse() const {
			LinearTransform<double> r;
			r._a = Real(1.0) / _a;
			r._b = -_b / _a;
			return r;
		}
	private:
		Real _a;
		Real _b;
	};

}