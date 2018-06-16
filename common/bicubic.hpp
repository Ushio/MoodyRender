#pragma once

#include <functional>
#include <glm/glm.hpp>
#include <glm/ext.hpp>

inline double bicubic_kernel(double t, double tt, double ttt, double f0, double f1, double f2, double f3) {
	glm::dmat4 m = {
		0.0, -1.0, 2.0, -1.0,
		2.0, 0.0, -5.0, 3.0,
		0.0, 1.0, 4.0, -3.0,
		0.0, 0.0, -1.0, 1.0
	};
	return 0.5 * glm::dot(glm::dvec4(1.0, t, tt, ttt), m * glm::dvec4(f0, f1, f2, f3));
}
inline double bicubic_kernel(double t, double f0, double f1, double f2, double f3) {
	double tt = t * t;
	double ttt = tt * t;
	return bicubic_kernel(t, tt, ttt, f0, f1, f2, f3);
}

// f(0.0) = sample(0)
// f(1.0) = sample(size - 1)
inline double bicubic_1d(double x /* 0.0 => 1.0 */, int size, std::function<double(int)> sample) {
	double index_f = x * (size - 1);
	int index1 = (int)std::floor(index_f);

	double values[4];
	for (int i = 0; i < 4; ++i) {
		int index = index1 - 1 + i;
		index = glm::clamp(index, 0, size - 1);
		values[i] = sample(index);
	}

	double t = index_f - index1;
	return bicubic_kernel(t, values[0], values[1], values[2], values[3]);
}

inline double bicubic_2d(double x /* 0.0 => 1.0 */, double y /* 0.0 => 1.0 */, int sizex, int sizey, std::function<double(int, int)> sample) {
	double index_xf = x * (sizex - 1);
	double index_yf = y * (sizey - 1);
	int index1x = (int)std::floor(index_xf);
	int index1y = (int)std::floor(index_yf);
	double tx = index_xf - index1x;
	double ty = index_yf - index1y;

	double values[4][4];
	for (int y = 0; y < 4; ++y) {
		int yi = index1y - 1 + y;
		yi = glm::clamp(yi, 0, sizey - 1);
		for (int x = 0; x < 4; ++x) {
			int xi = index1x - 1 + x;
			xi = glm::clamp(xi, 0, sizex - 1);
			values[y][x] = sample(xi, yi);
		}
	}

	double ttx = tx * tx;
	double tttx = ttx * tx;
	return bicubic_kernel(ty,
		bicubic_kernel(tx, ttx, tttx, values[0][0], values[0][1], values[0][2], values[0][3]),
		bicubic_kernel(tx, ttx, tttx, values[1][0], values[1][1], values[1][2], values[1][3]),
		bicubic_kernel(tx, ttx, tttx, values[2][0], values[2][1], values[2][2], values[2][3]),
		bicubic_kernel(tx, ttx, tttx, values[3][0], values[3][1], values[3][2], values[3][3])
	);
}