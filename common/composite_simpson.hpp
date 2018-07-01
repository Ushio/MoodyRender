#pragma once

#include <functional>

namespace rt {
	/*
	合成シンプソン公式
	nは偶数
	*/
	template <class Real, class F>
	inline Real composite_simpson(F f, int n, Real a, Real b) {
		assert(n % 2 == 0 && "n is must be even");

		Real sum = 0;
		Real h = (b - a) / n;
		for (int i = 1; i < n; ++i) {
			Real c = (i & 0x1) ? Real(4.0) : Real(2.0);
			Real x = a + h * i;
			sum += c * f(x);
		}
		sum += f(a) + f(b);
		return sum * h * Real(1.0 / 3.0);
	}
}