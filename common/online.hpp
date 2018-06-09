#pragma once
#include <vector>

namespace rt {
	/*
	https://qpp.bitbucket.io/post/variance/
	*/
	template <class Real>
	class OnlineMean {
	public:
		void addSample(Real x) {
			_mean = (x - _mean) / Real(_n + Real(1.0)) + _mean;
			_n++;
		}
		Real mean() const {
			return _mean;
		}
		int sampleCount() const {
			return _n;
		}
	private:
		int _n = 0;
		Real _mean = Real(0.0);
	};

	template <class Real>
	class OnlineVariance {
	public:
		void addSample(Real x) {
			Real mu_pre = _mean.mean();
			_mean.addSample(x);
			Real mu_new = _mean.mean();
			_m += (x - mu_pre) * (x - mu_new);
		}
		Real mean() const {
			return _mean.mean();
		}
		Real variance() const {
			return _mean.sampleCount() == 0 ? Real(0.0) : _m / _mean.sampleCount();
		}
		int sampleCount() const {
			return _mean.sampleCount();
		}
	private:
		OnlineMean<Real> _mean;
		Real _m = Real(0.0);
	};

	// Reference
	template <class Real>
	inline void mean_and_variance(const std::vector<Real> &xs, Real *mean, Real *variance) {
		auto sum = std::accumulate(xs.begin(), xs.end(), Real(0.0), std::plus<Real>());
		*mean = sum / xs.size();

		Real vsum = Real(0.0);
		for (Real x : xs) {
			Real s = *mean - x;
			vsum += s * s;
		}

		*variance = vsum / xs.size();
	}
}