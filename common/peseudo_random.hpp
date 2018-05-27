#pragma once
#include <algorithm>

// 基本的な擬似乱数
namespace rt {
	struct PeseudoRandom {
		virtual ~PeseudoRandom() {}

		virtual uint32_t generate() = 0;
		virtual double uniform() = 0;
		virtual double uniform(double a, double b) = 0;

		virtual float uniformf() { return (float)uniform(); }
		virtual float uniformf(float a, float b) { return (float)(uniform(a, b)); }
	};
	struct Xor : public PeseudoRandom {
		Xor() {

		}
		Xor(uint32_t seed) {
			_y = std::max(seed, 1u);
		}

		// 0 <= x <= 0x7FFFFFFF
		uint32_t generate() {
			_y = _y ^ (_y << 13); _y = _y ^ (_y >> 17);
			uint32_t value = _y = _y ^ (_y << 5); // 1 ~ 0xFFFFFFFF(4294967295
			return value >> 1;
		}
		// 0.0 <= x < 1.0
		double uniform() {
			return double(generate()) / double(0x80000000);
		}
		double uniform(double a, double b) {
			return a + (b - a) * double(uniform());
		}
	public:
		uint32_t _y = 2463534242;
	};
}