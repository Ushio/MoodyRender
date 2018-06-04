#pragma once

#include <glm/glm.hpp>
#include <glm/ext.hpp>
#include "coordinate.hpp"

namespace rt {
	inline float chi_plus(float x) {
		return x < 0.0f ? 0.0f : 1.0f;
	}

	inline float D_Beckmann(const glm::vec3 &n, const glm::vec3 &h, float alpha) {
		float cosTheta = glm::dot(n, h);
		if (cosTheta < 1.0e-5) {
			return 0.0f;
		}
		float cosTheta2 = cosTheta * cosTheta;
		float cosTheta4 = cosTheta2 * cosTheta2;
		float alpha2 = alpha * alpha;
		float chi = chi_plus(cosTheta);

		// \tan { \theta  } =\pm \frac { \sqrt { 1-\cos { \theta  }  }  }{ \cos { \theta  }  } \\ \tan ^{ 2 }{ \theta  } =\frac { 1-\cos ^{ 2 }{ \theta  }  }{ \cos ^{ 2 }{ \theta  }  } 
		float tanTheta2 = (1.0f - cosTheta2) / cosTheta2;
		return chi * std::exp(-tanTheta2 / alpha2) / (glm::pi<float>() * alpha2 * cosTheta4);
	}
	inline float lambda_beckmann(float cosTheta, float alpha) {
		float tanThetaO = std::sqrt(1.0f - cosTheta * cosTheta) / cosTheta;
		float a = 1.0f / (alpha * tanThetaO);
		return (std::erf(a) - 1.0f) * 0.5f + std::exp(-a * a) / (2.0f * a * std::sqrt(glm::pi<float>()));
	}
	inline float G2_height_correlated_beckmann(const glm::vec3 &omega_i, const glm::vec3 &omega_o, const glm::vec3 &omega_h, const glm::vec3 &n, float alpha) {
		float numer = chi_plus(glm::dot(omega_o, omega_h)) * chi_plus(glm::dot(omega_i, omega_h));
		float denom = (1.0f + lambda_beckmann(glm::dot(omega_o, n), alpha) + lambda_beckmann(glm::dot(omega_i, n), alpha));
		return numer / denom;
	}

	inline float G2_v_cavity(glm::vec3 L, glm::vec3 V, glm::vec3 H, glm::vec3 N) {
		float a = 2.0f * glm::dot(N, H) * glm::dot(N, V) / glm::dot(V, H);
		float b = 2.0f * glm::dot(N, H) * glm::dot(N, L) / glm::dot(L, H);
		return glm::min(glm::min(a, b), 1.0f);
	}

	//inline double G_kalemen(Vec3 L, Vec3 V, Vec3 H, Vec3 N, double alpha) {
	//	return 2.0 * glm::dot(N, L) * glm::dot(N, V) / (1.0 + glm::dot(L, V));
	//}

	inline float BeckmannMicrofacetBRDF_without_F(const glm::vec3 &omega_i, const glm::vec3 &omega_o, const glm::vec3 &omega_h, const glm::vec3 &Ng, float alpha) {
		float d = D_Beckmann(Ng, omega_h, alpha);
		float g = G2_height_correlated_beckmann(omega_i, omega_o, omega_h, Ng, alpha);

		float cos_term_wo = glm::dot(Ng, omega_o);
		float cos_term_wi = glm::dot(Ng, omega_i);

		return d * g / (4.0f * cos_term_wo * cos_term_wi);
	}

	struct NDFImportanceSampler {
		static glm::vec3 sample_wi_Beckmann(PeseudoRandom *random, float alpha, glm::vec3 wo, glm::vec3 Ng) {
			float theta = std::atan(std::sqrt(-alpha * alpha * std::log(1.0f - random->uniform())));
			float phi = random->uniform(0.0f, glm::two_pi<double>());
			glm::vec3 sample = polar_to_cartesian(theta, phi);

			glm::vec3 harf = from_bxdf(Ng, sample);
			glm::vec3 wi = glm::reflect(-wo, harf);
			return wi;
		}
		static float pdfBeckmann(glm::vec3 sampled_wi, float alpha, glm::vec3 wo, glm::vec3 Ng) {
			if (glm::dot(Ng, sampled_wi) <= 0.0f) {
				return std::numeric_limits<float>::max();
			}
			glm::vec3 half = glm::normalize(sampled_wi + wo);
			return D_Beckmann(Ng, half, alpha) * glm::dot(Ng, half) / (4.0f * glm::dot(sampled_wi, half));
		}
	};

	inline float fresnel_v(float n, float k, float cosTheta) {
		float n2_add_k2 = n * n + k * k;
		float numer = n2_add_k2 - 2.0f * n * cosTheta + cosTheta * cosTheta;
		float denom = n2_add_k2 + 2.0f * n * cosTheta + cosTheta * cosTheta;
		return numer / denom;
	}

	inline float fresnel_h(float n, float k, float cosTheta) {
		float n2_add_k2_cosTheta2 = (n * n + k * k) * cosTheta * cosTheta;
		float numer = n2_add_k2_cosTheta2 - 2.0f * n * cosTheta + 1.0f;
		float denom = n2_add_k2_cosTheta2 + 2.0f * n * cosTheta + 1.0f;
		return numer / denom;
	}

	inline float fresnel_unpolarized(float n, float k, float cosTheta) {
		return (fresnel_v(n, k, cosTheta) + fresnel_h(n, k, cosTheta)) * 0.5f;
	}

	// Air to Glass
	inline float fresnel_dielectrics(float cosTheta) {
		auto sqr = [](float x) { return x * x; };

		float eta_t = 1.5f; // for Glass
		float eta_i = 1.0f;
		float c = cosTheta;
		float g = std::sqrt(eta_t * eta_t / sqr(eta_i) - 1.0f + sqr(c));

		float a = 0.5f * sqr(g - c) / sqr(g + c);
		float b = 1.0f + sqr(c * (g + c) - 1.0f) / sqr(c * (g - c) + 1.0f);
		return a * b;
	}

	float fresnel_shlick(float f0, float cosTheta) {
		return f0 + (1.0f - f0) * std::pow(1.0f - cosTheta, 5);
	}


}

