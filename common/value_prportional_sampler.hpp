#pragma once

#include <vector>

namespace rt{
	template <class Real>
	class ValueProportionalSampler {
	public:
		ValueProportionalSampler() {}
		ValueProportionalSampler(const std::vector<Real> &values) {
			Real sumValue = Real(0.0);
			for (int i = 0; i < values.size(); ++i) {
				sumValue += values[i];
				_cumulativeAreas.push_back(sumValue);
			}
			_sumValue = sumValue;
			_values = values;
		}
		template <class T, class Func>
		ValueProportionalSampler(const std::vector<T> &values, Func valueAccess) {
			Real sumValue = Real(0.0);
			for (int i = 0; i < values.size(); ++i) {
				Real value = valueAccess(values[i]);
				sumValue += value;
				_values.push_back(value);
				_cumulativeAreas.push_back(sumValue);
			}
			_sumValue = sumValue;
		}

		int sample(PeseudoRandom *random) const {
			Real area_at = (Real)random->uniform(0.0, _sumValue);
			auto it = std::upper_bound(_cumulativeAreas.begin(), _cumulativeAreas.end(), area_at);
			std::size_t index = std::distance(_cumulativeAreas.begin(), it);
			index = std::min(index, _cumulativeAreas.size() - 1);
			return (int)index;
		}
		Real sumValue() const {
			return _sumValue;
		}
		Real probability(int index) const {
			return _values[index] / _sumValue;
		}
		int size() const {
			return (int)_values.size();
		}
	private:
		Real _sumValue = Real(0.0);
		std::vector<Real> _values;
		std::vector<Real> _cumulativeAreas;
	};
}