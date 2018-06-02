#pragma once

#include <functional>
#include <tuple>

namespace rt {
	// Simpson's Rule with cache
	template <class Real>
	struct SimpsonRange {
		SimpsonRange() {}
		SimpsonRange(std::function<Real(Real)> f, Real a, Real b) :
			_f(f),
			_a(a),
			_fa(f(a)),
			_m((a + b) * Real(0.5)),
			_fm(f(_m)),
			_b(b),
			_fb(f(b))
		{

		}
		SimpsonRange(std::function<Real(Real)> f, Real a, Real fa, Real b, Real fb) :
			_f(f),
			_a(a),
			_fa(fa),
			_m((a + b) * Real(0.5)),
			_fm(f(_m)),
			_b(b),
			_fb(fb)
		{
			_integral = integrate();
		}
		std::tuple<SimpsonRange, SimpsonRange> separate() const {
			return std::make_tuple(
				SimpsonRange(_f, _a, _fa, _m, _fm),
				SimpsonRange(_f, _m, _fm, _b, _fb)
			);
		}

		Real a() const {
			return _a;
		}
		Real fa() const {
			return _fa;
		}
		Real m() const {
			return _m;
		}
		Real fm() const {
			return _fm;
		}
		Real b() const {
			return _b;
		}
		Real fb() const {
			return _fb;
		}
		Real integral() const {
			return _integral;
		}

		Real evaluateApproximation(Real x) const {
			Real h = std::abs(_a - _b) * Real(0.5);
			Real a = (_fa + _fb - Real(2.0) * _fm) / (Real(2.0) * h * h);
			Real b = (_fb - _fa) / (Real(2.0) * h);
			Real c = _fm;
			auto sqr = [](Real x) { return x * x; };
			return a * sqr(x - _m) + b * (x - _m) + c;
		}
	private:
		inline Real integrate() {
			return (_b - _a) * Real(1.0 / 6.0) * (_fa + Real(4.0) * _fm + _fb);
		}
		std::function<Real(Real)> _f;
		Real _a;
		Real _fa;
		Real _m;
		Real _fm;
		Real _b;
		Real _fb;
		Real _integral;
	};

	template <class Real>
	Real adaptive_simpson(const SimpsonRange<Real> &simpsonRange, Real eps) {
		SimpsonRange<Real> L;
		SimpsonRange<Real> R;
		std::tie(L, R) = simpsonRange.separate();
		auto L_R_minus_integrate_a_b = L.integral() + R.integral() - simpsonRange.integral();
		if (std::abs(L_R_minus_integrate_a_b) <= Real(15.0) * eps) {
			return L.integral() + R.integral() + L_R_minus_integrate_a_b * Real(1.0 / 15.0);
		}
		return adaptive_simpson(L, eps) + adaptive_simpson(R, eps);
	}
}
