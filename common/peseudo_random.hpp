#pragma once
#include <algorithm>
#include <random>
#include <glm/glm.hpp>
#include <glm/ext.hpp>

// 基本的な擬似乱数
namespace rt {
	struct PeseudoRandom {
		virtual ~PeseudoRandom() {}

		/* double */

		// 0.0 <= x < 1.0
		virtual double uniform64f() = 0;

		// 0.0 <= x < 1.0
		double uniform() {
			return uniform64f();
		}
		// a <= x < b
		double uniform(double a, double b) {
			return glm::mix(a, b, uniform64f());
		}

		/* float */
		// 0.0 <= x < 1.0
		virtual float uniform32f() {
			return uniform64f();
		}

		// 0.0 <= x < 1.0
		float uniformf() {
			return uniform32f();
		}
		// a <= x < b
		float uniformf(float a, float b) {
			return glm::mix(a, b, uniform32f());
		}
	};

	// copy and paste from:
	//     https://ja.wikipedia.org/wiki/Xorshift
	struct Xor64 : public PeseudoRandom {
		Xor64() {

		}
		Xor64(uint64_t seed) {
			_x = std::max(seed, 1ULL);
		}
		uint64_t next() {
			_x = _x ^ (_x << 13);
			_x = _x ^ (_x >> 7);
			_x = _x ^ (_x << 17);
			return _x;
		}

		// copy and paste from:
		//     http://xoshiro.di.unimi.it/
		// より直観的な説明
		//     http://marupeke296.com/TIPS_No16_flaotrandom.html
		double uniform64f() override {
			uint64_t x = next();
			uint64_t bits = (0x3FFULL << 52) | (x >> 12);
			double value = *reinterpret_cast<double *>(&bits) - 1.0;
			return value;
		}
		float uniform32f() override {
			uint64_t x = next();
			uint32_t bits = ((uint32_t)x >> 9) | 0x3f800000;
			float value = *reinterpret_cast<float *>(&bits) - 1.0f;
			return value;
		}

		uint64_t _x = 88172645463325252ULL;
	};

	/*
	http://xoshiro.di.unimi.it/xoroshiro128plus.c
	*/
	struct XoroshiroPlus128 : public PeseudoRandom {
		XoroshiroPlus128() {
			splitmix sp;
			sp.x = 38927482;
			s[0] = std::max(sp.next(), 1ULL);
			s[1] = std::max(sp.next(), 1ULL);
		}
		XoroshiroPlus128(uint64_t seed) {
			splitmix sp;
			sp.x = seed;
			s[0] = std::max(sp.next(), 1ULL);
			s[1] = std::max(sp.next(), 1ULL);
		}

		double uniform64f() override {
			uint64_t x = next();
			uint64_t bits = (0x3FFULL << 52) | (x >> 12);
			return *reinterpret_cast<double *>(&bits) - 1.0;
		}
		float uniform32f() override {
			uint64_t x = next();
			uint32_t bits = ((uint32_t)x >> 9) | 0x3f800000;
			float value = *reinterpret_cast<float *>(&bits) - 1.0f;
			return value;
		}
		/* This is the jump function for the generator. It is equivalent
		to 2^64 calls to next(); it can be used to generate 2^64
		non-overlapping subsequences for parallel computations. */
		void jump() {
			static const uint64_t JUMP[] = { 0xdf900294d8f554a5, 0x170865df4b3201fc };

			uint64_t s0 = 0;
			uint64_t s1 = 0;
			for (int i = 0; i < sizeof JUMP / sizeof *JUMP; i++)
			{
				for (int b = 0; b < 64; b++) {
					if (JUMP[i] & UINT64_C(1) << b) {
						s0 ^= s[0];
						s1 ^= s[1];
					}
					next();
				}
			}

			s[0] = s0;
			s[1] = s1;
		}
	private:
		// http://xoshiro.di.unimi.it/splitmix64.c
		// for generate seed
		struct splitmix {
			uint64_t x = 0; /* The state can be seeded with any value. */
			uint64_t next() {
				uint64_t z = (x += 0x9e3779b97f4a7c15);
				z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9;
				z = (z ^ (z >> 27)) * 0x94d049bb133111eb;
				return z ^ (z >> 31);
			}
		};
		uint64_t rotl(const uint64_t x, int k) const {
			return (x << k) | (x >> (64 - k));
		}
		uint64_t next() {
			const uint64_t s0 = s[0];
			uint64_t s1 = s[1];
			const uint64_t result = s0 + s1;

			s1 ^= s0;
			s[0] = rotl(s0, 24) ^ s1 ^ (s1 << 16); // a, b
			s[1] = rotl(s1, 37); // c

			return result;
		}
	private:
		uint64_t s[2];
	};

	struct MT : public PeseudoRandom {
		MT() {

		}
		MT(uint64_t seed) :_engine(seed) {

		}
		double uniform64f() override {
			std::uniform_real_distribution<> d(0.0, 1.0);
			return d(_engine);
		}
		std::mt19937 _engine;
	};
}