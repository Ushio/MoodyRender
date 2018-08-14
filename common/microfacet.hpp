#pragma once

#include <glm/glm.hpp>
#include <glm/ext.hpp>

#include "peseudo_random.hpp"
#include "coordinate.hpp"
#include "composite_simpson.hpp"
#include "serializable_buffer.hpp"
#include "value_prportional_sampler.hpp"
#include "composite_simpson.hpp"

namespace rt {
	inline double chi_plus(double x) {
		return x <= 0.0 ? 0.0 : 1.0;
	}

	inline double D_Beckmann(const glm::dvec3 &n, const glm::dvec3 &h, double alpha) {
		double cosTheta = glm::dot(n, h);

		// chi+
		if (cosTheta < 1.0e-5) {
			return 0.0;
		}

		double cosTheta2 = cosTheta * cosTheta;
		double cosTheta4 = cosTheta2 * cosTheta2;
		double alpha2 = alpha * alpha;
		double chi = chi_plus(cosTheta);

		// \tan { \theta  } =\pm \frac { \sqrt { 1-\cos { \theta  }  }  }{ \cos { \theta  }  } \\ \tan ^{ 2 }{ \theta  } =\frac { 1-\cos ^{ 2 }{ \theta  }  }{ \cos ^{ 2 }{ \theta  }  } 
		double tanTheta2 = (1.0 - cosTheta2) / cosTheta2;
		return chi * std::exp(-tanTheta2 / alpha2) / (glm::pi<double>() * alpha2 * cosTheta4);
	}
	inline double lambda_beckmann(double cosTheta, double alpha) {
		double tanThetaO = std::sqrt(1.0 - cosTheta * cosTheta) / cosTheta;
		double a = 1.0 / (alpha * tanThetaO);
		return (std::erf(a) - 1.0) * 0.5 + std::exp(-a * a) / (2.0 * a * std::sqrt(glm::pi<double>()));
	}
	inline double G2_height_correlated_beckmann(const glm::dvec3 &omega_i, const glm::dvec3 &omega_o, const glm::dvec3 &omega_h, const glm::dvec3 &n, double alpha) {
		double numer = chi_plus(glm::dot(omega_o, omega_h)) * chi_plus(glm::dot(omega_i, omega_h));
		double denom = (1.0 + lambda_beckmann(glm::dot(omega_o, n), alpha) + lambda_beckmann(glm::dot(omega_i, n), alpha));
		return numer / denom;
	}

	inline double G2_v_cavity(glm::dvec3 L, glm::dvec3 V, glm::dvec3 H, glm::dvec3 N) {
		double a = 2.0 * glm::dot(N, H) * glm::dot(N, V) / glm::max(glm::dot(V, H), 0.0);
		double b = 2.0 * glm::dot(N, H) * glm::dot(N, L) / glm::max(glm::dot(L, H), 0.0);
		return glm::min(glm::min(a, b), 1.0);
	}

	inline double G1_v_cavity(glm::dvec3 omega, glm::dvec3 H, glm::dvec3 N) {
		double a = 2.0 * glm::dot(N, H) * glm::dot(N, omega) / glm::max(glm::dot(omega, H), 0.0);
		return glm::min(a, 1.0);
	}

	//inline double G_kalemen(Vec3 L, Vec3 V, Vec3 H, Vec3 N, double alpha) {
	//	return 2.0 * glm::dot(N, L) * glm::dot(N, V) / (1.0 + glm::dot(L, V));
	//}

	//inline double beckmannMicrofacetBRDF_without_F(const glm::dvec3 &omega_i, const glm::dvec3 &omega_o, const glm::dvec3 &omega_h, const glm::dvec3 &Ng, double alpha) {
	//	double d = D_Beckmann(Ng, omega_h, alpha);
	//	double g = G2_height_correlated_beckmann(omega_i, omega_o, omega_h, Ng, alpha);

	//	double cos_term_wo = glm::dot(Ng, omega_o);
	//	double cos_term_wi = glm::dot(Ng, omega_i);

	//	return chi_plus(glm::dot(Ng, omega_i)) * chi_plus(glm::dot(Ng, omega_o)) *  d * g / (4.0 * cos_term_wo * cos_term_wi);
	//}
	
	struct BeckmannMicrosurfaceImportanceSampler {
		static glm::dvec3 sample(PeseudoRandom *random, double alpha) {
			double theta = std::atan(std::sqrt(-alpha * alpha * std::log(random->uniform())));
			double phi = random->uniform(0.0, glm::two_pi<double>());
			return polar_to_cartesian(theta, phi);
		}
		static double pdf(glm::dvec3 sampled_wi, double alpha, glm::dvec3 wo, glm::dvec3 Ng, glm::dvec3 *m = nullptr) {
			glm::dvec3 half = glm::normalize(sampled_wi + wo);
			if (m) {
				*m = half;
			}
			return D_Beckmann(Ng, half, alpha) * glm::dot(Ng, half);
		}
	};

	struct BeckmannImportanceSampler {
		// サンプリング範囲が半球ではないことに注意
		static glm::dvec3 sample(PeseudoRandom *random, double alpha, glm::dvec3 wo, glm::dvec3 Ng) {
			glm::dvec3 sample = BeckmannMicrosurfaceImportanceSampler::sample(random, alpha);

			ArbitraryBRDFSpace basis(Ng);
			glm::dvec3 h = basis.localToGlobal(sample);
			glm::dvec3 wi = glm::reflect(-wo, h);
			return wi;
		}
		// 裏面はサポートしない
		static double pdf(glm::dvec3 sampled_wi, double alpha, glm::dvec3 wo, glm::dvec3 Ng) {
			glm::dvec3 m;
			double pdf_m = BeckmannMicrosurfaceImportanceSampler::pdf(sampled_wi, alpha, wo, Ng, &m);

			// glm::dot(sampled_wi, half)が0になるのは、
			// wiとwoが正反対の向き、つまりかならず裏側であるので、普段は問題にならない
			return pdf_m / (4.0 * glm::dot(sampled_wi, m));
		}
	};

	struct VCavityBeckmannVisibleNormalSampler {
		static glm::dvec3 sample(PeseudoRandom *random, double alpha, glm::dvec3 wo, glm::dvec3 Ng) {
			double phi = random->uniform(0.0, glm::two_pi<double>());

			glm::dvec3 omega_m = BeckmannMicrosurfaceImportanceSampler::sample(random, alpha);
			glm::dvec3 omega_m_dot = glm::dvec3(-omega_m.x, -omega_m.y, omega_m.z);

			ArbitraryBRDFSpace basis(Ng);
			glm::dvec3 wo_local = basis.globalToLocal(wo);

			double visible     = glm::max(glm::dot(wo_local, omega_m),     0.0);
			double visible_dot = glm::max(glm::dot(wo_local, omega_m_dot), 0.0);
			double u = visible_dot / (visible + visible_dot);

			glm::dvec3 sample = random->uniform() < u ? omega_m_dot : omega_m;

			glm::dvec3 h = basis.localToGlobal(sample);
			glm::dvec3 wi = glm::reflect(-wo, h);
			return wi;
		}

		// 裏面はサポートしない
		static double pdf(glm::dvec3 sampled_wi, double alpha, glm::dvec3 wo, glm::dvec3 Ng) {
			glm::dvec3 wm = glm::normalize(sampled_wi + wo);
			double cosThetaO = glm::dot(wo, Ng);
			return G1_v_cavity(wo, wm, Ng) * glm::max(glm::dot(wo, wm), 0.0) * D_Beckmann(Ng, wm, alpha) / (cosThetaO * (4.0 * glm::dot(sampled_wi, wm)));
		}
	};

	inline double velvet_D(const glm::dvec3 &n, const glm::dvec3 &h, double r) {
		double cosTheta = glm::dot(n, h);
		if (cosTheta < 0.0) {
			return 0.0;
		}
		// double sinTheta = std::sin(std::acos(cosTheta));
		double sinTheta = std::sqrt(std::max(1.0 - cosTheta * cosTheta, 0.0));
		return (2.0 + 1.0 / r) * std::pow(sinTheta, 1.0 / r) / (glm::pi<double>() * 2.0);
	}

	// a, b, c, d, e
	// power_of_one_minus_r = (1-r)^2
	inline double velvet_params_interpolate(int i, double power_of_one_minus_r) {
		static const double p0[5] = { 25.3245, 3.32435, 0.16801, -1.27393, -4.85967 };
		static const double p1[5] = { 21.5473, 3.82987, 0.19823, -1.97760, -4.32054 };
		// k * p0 + (1 - k) * p1
		// k * p0 + p1 - k * p1
		// p1 + k * (p0 - p1);
		return glm::mix(p1[i], p0[i], power_of_one_minus_r);
	}
	inline double velvet_L(double x, double r) {
		double one_minus_r = 1.0 - r;
		double power_of_one_minus_r = one_minus_r * one_minus_r;
		double a = velvet_params_interpolate(0, power_of_one_minus_r);
		double b = velvet_params_interpolate(1, power_of_one_minus_r);
		double c = velvet_params_interpolate(2, power_of_one_minus_r);
		double d = velvet_params_interpolate(3, power_of_one_minus_r);
		double e = velvet_params_interpolate(4, power_of_one_minus_r);
		return a / (1.0 + b * std::pow(x, c)) + d * x + e;
	}
	inline double velvet_lambda(double cosTheta, double r) {
		if (cosTheta < 0.5) {
			return std::exp(velvet_L(cosTheta, r));
		}
		return std::exp(2.0 * velvet_L(0.5, r) - velvet_L(1.0 - cosTheta, r));
	}
	//inline double velvet_lambda_dot(double cosTheta, double r) {
	//	return std::pow(velvet_lambda(cosTheta, r), 1.0 + 2.0 * std::pow(1.0 - cosTheta, 8));
	//}
	inline double velvet_G1(double cosTheta, double r) {
		return chi_plus(cosTheta) / (1.0 + velvet_lambda(cosTheta, r));
	}
	inline double velvet_G2(double cosThetaO, double cosThetaI, double r) {
		return chi_plus(cosThetaO) * chi_plus(cosThetaI) / (1.0 + velvet_lambda(cosThetaO, r) + velvet_lambda(cosThetaI, r));
	}
	//inline double velvet_G1_dot(double cosTheta, double r) {
	//	return chi_plus(cosTheta) / (1.0 + velvet_lambda_dot(cosTheta, r));
	//}
	//inline double velvet_G2_dot(double cosThetaO, double cosThetaI, double r) {
	//	return chi_plus(cosThetaO) * chi_plus(cosThetaI) / (1.0 + velvet_lambda_dot(cosThetaO, r) + velvet_lambda_dot(cosThetaI, r));
	//}
	// まだうまく動作しない
	struct VelvetSampler {
		// サンプリング範囲が半球ではないことに注意
		static glm::dvec3 sample(PeseudoRandom *random, double alpha, glm::dvec3 wo, glm::dvec3 Ng) {
			double phi = random->uniform(0.0, glm::two_pi<double>());

			glm::dvec3 omega_m;
			{
				double u = random->uniform();
				double theta = std::asin(std::pow(u, alpha / (alpha * 2.0 + 1.0)));
				double phi = random->uniform(0.0, glm::two_pi<double>());
				omega_m = polar_to_cartesian(theta, phi);
			}
			//ArbitraryBRDFSpace basis(Ng);
			//return basis.localToGlobal(omega_m);
			ArbitraryBRDFSpace basis(Ng);
			glm::dvec3 h = basis.localToGlobal(omega_m);
			glm::dvec3 wi = glm::reflect(-wo, h);
			return wi;
		}
		// 裏面はサポートしない
		static double pdf(glm::dvec3 sampled_wi, double alpha, glm::dvec3 wo, glm::dvec3 Ng) {
			glm::dvec3 half = glm::normalize(sampled_wi + wo);
			double pdf_m = velvet_D(Ng, half, alpha) * glm::dot(Ng, half);

			// glm::dot(sampled_wi, half)が0になるのは、
			// wiとwoが正反対の向き、つまりかならず裏側であるので、普段は問題にならない
			return pdf_m / (4.0 * glm::dot(sampled_wi, half));
		}
	};

	// fresnel conductor
	inline double fresnel_v(double n, double k, double cosTheta) {
		double n2_add_k2 = n * n + k * k;
		double numer = n2_add_k2 - 2.0 * n * cosTheta + cosTheta * cosTheta;
		double denom = n2_add_k2 + 2.0 * n * cosTheta + cosTheta * cosTheta;
		return numer / denom;
	}

	inline double fresnel_h(double n, double k, double cosTheta) {
		double n2_add_k2_cosTheta2 = (n * n + k * k) * cosTheta * cosTheta;
		double numer = n2_add_k2_cosTheta2 - 2.0 * n * cosTheta + 1.0;
		double denom = n2_add_k2_cosTheta2 + 2.0 * n * cosTheta + 1.0;
		return numer / denom;
	}

	inline double fresnel_unpolarized(double n, double k, double cosTheta) {
		return (fresnel_v(n, k, cosTheta) + fresnel_h(n, k, cosTheta)) * 0.5;
	}

	// 1.5: grass
	inline double fresnel_dielectrics(double cosTheta, double eta_t, double eta_i) {
		auto sqr = [](double x) { return x * x; };

		double c = cosTheta;
		double g = std::sqrt(eta_t * eta_t / sqr(eta_i) - 1.0 + sqr(c));

		double a = 0.5 * sqr(g - c) / sqr(g + c);
		double b = 1.0 + sqr(c * (g + c) - 1.0) / sqr(c * (g - c) + 1.0);
		return a * b;
	}

	inline double fresnel_shlick(double f0, double cosTheta) {
		return f0 + (1.0 - f0) * std::pow(1.0 - cosTheta, 5);
	}

	// conductor fresnel avg
	inline double fresnel_avg(double n, double k) {
		//return rt::composite_simpson<double>([&](double theta) {
		//	double cosTheta = std::cos(theta);
		//	double sinTheta = std::sin(theta);
		//	return 2.0 * fresnel_unpolarized(n, k, cosTheta) * cosTheta * sinTheta;
		//}, 6, 0.0, glm::pi<double>() * 0.5);
		return 2.0 * rt::composite_simpson<double>([&](double cosTheta) {
			return fresnel_unpolarized(n, k, cosTheta) * cosTheta;
		}, 10, 0.0, 1.0);
	}

	inline double CoupledBRDF_Proportional(double theta, double alpha, std::function<double(double, double)> specularAlbedo) {
		double cosTheta = std::cos(theta);
		double sinTheta = std::sin(theta);
		return (1.0 - specularAlbedo(alpha, cosTheta)) * cosTheta * sinTheta;
	}

	class CoupledBRDFSampler {
	public:
		// specularAlbedo(alpha, cosTheta)
		void build(std::function<double(double, double)> specularAlbedo) {
			int kAlphaCount = 32;
			int kSampleBlockCount = 128;

			_discreteSamplers.resize(kAlphaCount);
			for (int i = 0; i < kAlphaCount; ++i) {
				double alpha = indexToAlpha(i, kAlphaCount);

				std::vector<double> values(kSampleBlockCount);
				for (int j = 0; j < kSampleBlockCount; ++j) {
					values[j] = CoupledBRDF_Proportional(indexToTheta(j, kSampleBlockCount), alpha, specularAlbedo);
				}
				_discreteSamplers[i] = ValueProportionalSampler<double>(values);
			}
		}
		double sampleTheta(double alpha, PeseudoRandom *random) const {
			int alphaIndex = alphaToIndex(alpha, (int)_discreteSamplers.size());
			const ValueProportionalSampler<double> &sampler = _discreteSamplers[alphaIndex];
			int indexTheta = sampler.sample(random);
			auto thetaRange = indexToThetaRange(indexTheta, (int)sampler.size());
			return random->uniform(thetaRange.first, thetaRange.second);
		}
		double probability(double alpha, double theta) const {
			int alphaIndex = alphaToIndex(alpha, (int)_discreteSamplers.size());
			const ValueProportionalSampler<double> &sampler = _discreteSamplers[alphaIndex];
			int thetaIndex = thetaToIndex(theta, sampler.size());
			return sampler.probability(thetaIndex);
		}
		int thetaSize(double alpha) const {
			int alphaIndex = alphaToIndex(alpha, (int)_discreteSamplers.size());
			const ValueProportionalSampler<double> &sampler = _discreteSamplers[alphaIndex];
			return sampler.size();
		}
	private:
		// 0     0.5     1
		// |------|------|
		int alphaToIndex(double alpha, int n)const
		{
			int index = (int)(alpha * n);
			index = std::min(index, n - 1);
			index = std::max(index, 0);
			return index;
		}
		double indexToAlpha(int index, int n) const {
			double wide = 1.0 / n;
			return wide * 0.5 + index * wide;
		}

		// 0     0.5     1
		// |------|------|
		int thetaToIndex(double theta, int n) const
		{
			int index = (int)((theta * 2.0 / glm::pi<double>()) * n);
			index = std::min(index, n - 1);
			index = std::max(index, 0);
			return index;
		}
		double indexToTheta(int index, int n) const {
			double wide = glm::pi<double>() * 0.5 / n;
			return wide * 0.5 + index * wide;
		}
		std::pair<double, double> indexToThetaRange(int index, int n) const {
			double wide = glm::pi<double>() * 0.5 / n;
			return std::make_pair(index * wide, (index + 1) * wide);
		}

		// alpha => table
		std::vector<ValueProportionalSampler<double>> _discreteSamplers;
	};

	class CoupledBRDFConductor {
	public:
		static SpecularAlbedo & specularAlbedo() {
			static SpecularAlbedo s_specularAlbedo;
			return s_specularAlbedo;
		}
		static SpecularAvgAlbedo& specularAvgAlbedo() {
			static SpecularAvgAlbedo s_specularAvgAlbedo;
			return s_specularAvgAlbedo;
		}
		static CoupledBRDFSampler &sampler() {
			static CoupledBRDFSampler s_sampler;
			return s_sampler;
		}

		static void load(const char *specularAlbedoBin, const char *specularAvgAlbedoBin) {
			loadFromBinary(specularAlbedo(), specularAlbedoBin);
			loadFromBinary(specularAvgAlbedo(), specularAvgAlbedoBin);
			sampler().build([](double alpha, double cosTheta) {
				return rt::CoupledBRDFConductor::specularAlbedo().sample(alpha, cosTheta);
			});
		}
	};
	class CoupledBRDFDielectrics {
	public:
		static SpecularAlbedo & specularAlbedo() {
			static SpecularAlbedo s_specularAlbedo;
			return s_specularAlbedo;
		}
		static SpecularAvgAlbedo& specularAvgAlbedo() {
			static SpecularAvgAlbedo s_specularAvgAlbedo;
			return s_specularAvgAlbedo;
		}
		static CoupledBRDFSampler &sampler() {
			static CoupledBRDFSampler s_sampler;
			return s_sampler;
		}
		static void load(const char *specularAlbedoBin, const char *specularAvgAlbedoBin) {
			loadFromBinary(specularAlbedo(), specularAlbedoBin);
			loadFromBinary(specularAvgAlbedo(), specularAvgAlbedoBin);
			sampler().build([](double alpha, double cosTheta) {
				return rt::CoupledBRDFDielectrics::specularAlbedo().sample(alpha, cosTheta);
			});
		}
	};

	class CoupledBRDFVelvet {
	public:
		static SpecularAlbedo & specularAlbedo() {
			static SpecularAlbedo s_specularAlbedo;
			return s_specularAlbedo;
		}
		static SpecularAvgAlbedo& specularAvgAlbedo() {
			static SpecularAvgAlbedo s_specularAvgAlbedo;
			return s_specularAvgAlbedo;
		}
		static CoupledBRDFSampler &sampler() {
			static CoupledBRDFSampler s_sampler;
			return s_sampler;
		}
		static void load(const char *specularAlbedoBin, const char *specularAvgAlbedoBin) {
			loadFromBinary(specularAlbedo(), specularAlbedoBin);
			loadFromBinary(specularAvgAlbedo(), specularAvgAlbedoBin);
			sampler().build([](double alpha, double cosTheta) {
				return rt::CoupledBRDFVelvet::specularAlbedo().sample(alpha, cosTheta);
			});
		}
	};
}

