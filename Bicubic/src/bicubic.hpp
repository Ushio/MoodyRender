#pragma once

#include <functional>
#include <glm/glm.hpp>
#include <glm/ext.hpp>

inline float bicubic_kernel(float t, float tt, float ttt, float f0, float f1, float f2, float f3) {
	glm::mat4 m = {
		0.0f, -1.0f, 2.0f, -1.0f,
		2.0f, 0.0f, -5.0f, 3.0f,
		0.0f, 1.0f, 4.0f, -3.0f,
		0.0f, 0.0f, -1.0f, 1.0f
	};
	return 0.5f * glm::dot(glm::vec4(1.0f, t, tt, ttt), m * glm::vec4(f0, f1, f2, f3));
}
inline float bicubic_kernel(float t, float f0, float f1, float f2, float f3) {
	float tt = t * t;
	float ttt = tt * t;
	return bicubic_kernel(t, tt, ttt, f0, f1, f2, f3);
}

// f(0.0) = sample(0)
// f(1.0) = sample(size - 1)
inline float bicubic_1d(float x /* 0.0 => 1.0 */, int size, std::function<float(int)> sample) {
	float index_f = x * (size - 1);
	int index1 = (int)std::floor(index_f);

	float values[4];
	for (int i = 0; i < 4; ++i) {
		int index = index1 - 1 + i;
		index = glm::clamp(index, 0, size - 1);
		values[i] = sample(index);
	}

	float t = index_f - index1;
	return bicubic_kernel(t, values[0], values[1], values[2], values[3]);
}

inline float bicubic_2d(float x /* 0.0 => 1.0 */, float y /* 0.0 => 1.0 */, int sizex, int sizey, std::function<float(int, int)> sample) {
	float index_xf = x * (sizex - 1);
	float index_yf = y * (sizey - 1);
	int index1x = (int)std::floor(index_xf);
	int index1y = (int)std::floor(index_yf);
	float tx = index_xf - index1x;
	float ty = index_yf - index1y;

	float values[4][4];
	for (int y = 0; y < 4; ++y) {
		int yi = index1y - 1 + y;
		yi = glm::clamp(yi, 0, sizey - 1);
		for (int x = 0; x < 4; ++x) {
			int xi = index1x - 1 + x;
			xi = glm::clamp(xi, 0, sizex - 1);
			values[y][x] = sample(xi, yi);
		}
	}

	float ttx = tx * tx;
	float tttx = ttx * tx;
	return bicubic_kernel(ty,
		bicubic_kernel(tx, ttx, tttx, values[0][0], values[0][1], values[0][2], values[0][3]),
		bicubic_kernel(tx, ttx, tttx, values[1][0], values[1][1], values[1][2], values[1][3]),
		bicubic_kernel(tx, ttx, tttx, values[2][0], values[2][1], values[2][2], values[2][3]),
		bicubic_kernel(tx, ttx, tttx, values[3][0], values[3][1], values[3][2], values[3][3])
	);
}