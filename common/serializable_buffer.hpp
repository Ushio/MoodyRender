#pragma once

#include <vector>
#include <functional>
#include <fstream>

#include <cereal/cereal.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/archives/xml.hpp>

#include <tbb/tbb.h>

#include "bicubic.hpp"

namespace rt {
	class SpecularAlbedo {
	public:
		SpecularAlbedo() {}

		void save(const char *filename) {
			std::ofstream ofs(filename);
			{
				cereal::XMLOutputArchive o_archive(ofs);
				o_archive(*this);
			}
		}
		void load(const char *filename) {
			std::ifstream ifs(filename);
			{
				cereal::XMLInputArchive i_archive(ifs);
				i_archive(*this);
			}
		}

		// evaluate(alpha, cosTheta)
		void build(int alphaSize, int cosThetaSize, std::function<float(float, float)> evaluate) {
			_alphaSize = alphaSize;
			_cosThetaSize = cosThetaSize;
			_values.resize(_alphaSize * _cosThetaSize);

			tbb::parallel_for(tbb::blocked_range<int>(0, _cosThetaSize), [&](const tbb::blocked_range<int> &range) {
				for (int j = range.begin(); j < range.end(); ++j) {
					// cosTheta == 0 を回避するために j == 0 を回避する
					float cosTheta = (float)std::max(j, 1) / (_cosThetaSize - 1);
					for (int i = 0; i < _alphaSize; ++i) {
						// alpha == 0 を回避するために i == 0 を回避する
						float alpha = (float)std::max(i, 1) / (_alphaSize - 1);
						set(i, j, evaluate(alpha, cosTheta));
					}

					printf("%d line done.\n", j);
				}
			});
		}
		float sample(float alpha, float cosTheta) const {
			return bicubic_2d(alpha, cosTheta, _alphaSize, _cosThetaSize, [&](int x, int y) {
				return get(x, y);
			});
		}
		void set(int x, int y, float value) {
			_values[y * _alphaSize + x] = value;
		}
		float get(int x, int y) const {
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
		std::vector<float> _values;
	};

	class SpecularAvgAlbedo {
	public:
		SpecularAvgAlbedo() {}

		void save(const char *filename) {
			std::ofstream ofs(filename);
			{
				cereal::XMLOutputArchive o_archive(ofs);
				o_archive(*this);
			}
		}
		void load(const char *filename) {
			std::ifstream ifs(filename);
			{
				cereal::XMLInputArchive i_archive(ifs);
				i_archive(*this);
			}
		}

		// evaluate(alpha)
		void build(int alphaSize, std::function<float(float)> evaluate) {
			_alphaSize = alphaSize;
			_values.resize(_alphaSize);

			for (int i = 0; i < _alphaSize; ++i) {
				// alpha == 0 を回避するために i == 0 を回避する
				float alpha = (float)std::max(i, 1) / (_alphaSize - 1);
				set(i, evaluate(alpha));
			}
		}

		float sample(float alpha) const {
			return bicubic_1d(alpha, _alphaSize, [&](int x) {
				return get(x);
			});
		}
		void set(int x, float value) {
			_values[x] = value;
		}
		float get(int x) const {
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
		std::vector<float> _values;
	};


	class I_Dot_Inverse {
	public:
		I_Dot_Inverse() {}

		void save(const char *filename) {
			std::ofstream ofs(filename);
			{
				cereal::XMLOutputArchive o_archive(ofs);
				o_archive(*this);
			}
		}
		void load(const char *filename) {
			std::ifstream ifs(filename);
			{
				cereal::XMLInputArchive i_archive(ifs);
				i_archive(*this);
			}
		}

		// evaluate I_dot_inverse(alpha, u)
		void build(int alphaSize, int uSize, std::function<float(float, float)> I_dot_inverse) {
			_alphaSize = alphaSize;
			_uSize = uSize;
			_values.resize(_alphaSize * _uSize);

			tbb::parallel_for(tbb::blocked_range<int>(0, _uSize), [&](const tbb::blocked_range<int> &range) {
				for (int j = range.begin(); j < range.end(); ++j) {
					float u = (float)j / (_uSize - 1);
					for (int i = 0; i < _alphaSize; ++i) {
						// alpha == 0 を回避するために i == 0 を回避する
						float alpha = (float)std::max(i, 1) / (_alphaSize - 1);
						float value = I_dot_inverse(1.0 - alpha, u);
						set(i, j, value);
					}

					printf("%d line done\n", j);
				}
			});
		}
		float sample(float alpha, float u) const {
			return bicubic_2d(alpha, u, _alphaSize, _uSize, [&](int x, int y) {
				return get(x, y);
			});
		}
		void set(int x, int y, float value) {
			_values[y * _alphaSize + x] = value;
		}
		float get(int x, int y) const {
			return _values[y * _alphaSize + x];
		}
		int alphaSize() const {
			return _alphaSize;
		}
		int uSize() const {
			return _uSize;
		}
	private:
		friend class cereal::access;

		template<class Archive>
		void serialize(Archive & archive)
		{
			archive(CEREAL_NVP(_alphaSize), CEREAL_NVP(_uSize), CEREAL_NVP(_values));
		}

		int _alphaSize = 0;
		int _uSize = 0;
		std::vector<float> _values;
	};
}