#pragma once
#include <chrono>

// 
namespace rt {
	class Stopwatch {
	public:
		Stopwatch() :_beginAt(std::chrono::high_resolution_clock::now()) {
		}
		double elapsed() const {
			auto n = std::chrono::high_resolution_clock::now();
			return (double)std::chrono::duration_cast<std::chrono::microseconds>(n - _beginAt).count() * 0.001 * 0.001;
		}
	private:
		std::chrono::high_resolution_clock::time_point _beginAt;
	};
}