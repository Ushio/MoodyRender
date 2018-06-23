#pragma once

#include <vector>
#include <functional>
#include <fstream>

#include <cereal/cereal.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/archives/xml.hpp>
#include <cereal/archives/portable_binary.hpp>

#include <tbb/tbb.h>

#include "bicubic.hpp"

namespace rt {
	template <class T>
	inline void saveAsBinary(const T &value, const char *filename) {
		std::ofstream ofs(filename, std::ios::binary);
		{
			cereal::PortableBinaryOutputArchive o_archive(ofs);
			o_archive(value);
		}
	}
	template <class T>
	inline void loadFromBinary(T &value, const char *filename) {
		std::ifstream ifs(filename, std::ios::binary);
		{
			cereal::PortableBinaryInputArchive i_archive(ifs);
			i_archive(value);
		}
	}

	class SpecularAlbedo {
	public:
		SpecularAlbedo() {}

		// evaluate(alpha, cosTheta)
		void build(int alphaSize, int cosThetaSize, std::function<double(double, double)> evaluate) {
			_alphaSize = alphaSize;
			_cosThetaSize = cosThetaSize;
			_values.resize(_alphaSize * _cosThetaSize);

			tbb::parallel_for(tbb::blocked_range<int>(0, _cosThetaSize), [&](const tbb::blocked_range<int> &range) {
				for (int j = range.begin(); j < range.end(); ++j) {
					// cosTheta == 0 を回避するために j == 0 を回避する
					double cosTheta = (double)std::max(j, 1) / (_cosThetaSize - 1);
					for (int i = 0; i < _alphaSize; ++i) {
						// alpha == 0 を回避するために i == 0 を回避する
						double alpha = (double)std::max(i, 1) / (_alphaSize - 1);
						set(i, j, evaluate(alpha, cosTheta));
					}

					printf("%d line done.\n", j);
				}
			});
		}
		double sample(double alpha, double cosTheta) const {
			return bicubic_2d<double>(alpha, cosTheta, _alphaSize, _cosThetaSize, [&](int x, int y) {
				return get(x, y);
			});
		}
		void set(int x, int y, double value) {
			_values[y * _alphaSize + x] = value;
		}
		double get(int x, int y) const {
			return _values[y * _alphaSize + x];
		}
		int alphaSize() const {
			return _alphaSize;
		}
		int cosThetaSize() const {
			return _cosThetaSize;
		}
	private:
		friend class cereal::access;

		template<class Archive>
		void serialize(Archive & archive)
		{
			archive(CEREAL_NVP(_alphaSize), CEREAL_NVP(_cosThetaSize), CEREAL_NVP(_values));
		}

		int _alphaSize = 0;
		int _cosThetaSize = 0;
		std::vector<double> _values;
	};

	class SpecularAvgAlbedo {
	public:
		SpecularAvgAlbedo() {}

		// evaluate(alpha)
		void build(int alphaSize, std::function<double(double)> evaluate) {
			_alphaSize = alphaSize;
			_values.resize(_alphaSize);

			for (int i = 0; i < _alphaSize; ++i) {
				// alpha == 0 を回避するために i == 0 を回避する
				double alpha = (double)std::max(i, 1) / (_alphaSize - 1);
				set(i, evaluate(alpha));
			}
		}

		double sample(double alpha) const {
			return bicubic_1d<double>(alpha, _alphaSize, [&](int x) {
				return get(x);
			});
		}
		void set(int x, double value) {
			_values[x] = value;
		}
		double get(int x) const {
			return _values[x];
		}
		int alphaSize() const {
			return _alphaSize;
		}
	private:
		friend class cereal::access;

		template<class Archive>
		void serialize(Archive & archive)
		{
			archive(CEREAL_NVP(_alphaSize), CEREAL_NVP(_values));
		}
		int _alphaSize = 0;
		std::vector<double> _values;
	};
}