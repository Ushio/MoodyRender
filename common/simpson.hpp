#pragma once

namespace rt {
	/*
	合成シンプソン公式
	nは偶数
	*/
	template <int n, class F>
	inline double integrate_composite_simpson(const F &f, double a, double b) {
		static_assert(n % 2 == 0, "n is must be even");

		double sum = 0;
		double h = (b - a) / n;
		for (int i = 1; i < n; ++i) {
			double c = (i & 0x1) ? 4.0 : 2.0;
			double x = a + h * i;
			sum += c * f(x);
		}
		sum += f(a) + f(b);
		return sum * h / 3.0;
	}
}