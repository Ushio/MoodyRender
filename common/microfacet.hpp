﻿#pragma once

#include <glm/glm.hpp>
#include <glm/ext.hpp>

namespace rt {
	inline float chi_plus(float x) {
		return x < 0.0f ? 0.0f : 1.0f;
	}

	inline float D_Beckman(const glm::vec3 &n, const glm::vec3 &h, float alpha) {
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
	inline float G2_height_correlated_beckmann(const glm::vec3 &omega_i, const glm::vec3 &omega_o, const glm::vec3 &omega_h, const glm::vec3 &n, double alpha) {
		float numer = chi_plus(glm::dot(omega_o, omega_h)) * chi_plus(glm::dot(omega_i, omega_h));
		float denom = (1.0f + lambda_beckmann(glm::dot(omega_o, n), alpha) + lambda_beckmann(glm::dot(omega_i, n), alpha));
		return numer / denom;
	}

	inline float G_v_cavity(glm::vec3 L, glm::vec3 V, glm::vec3 H, glm::vec3 N, float alpha) {
		float a = 2.0f * glm::dot(N, H) * glm::dot(N, V) / glm::dot(V, H);
		float b = 2.0f * glm::dot(N, H) * glm::dot(N, L) / glm::dot(L, H);
		return glm::min(glm::min(a, b), 1.0f);
	}

	//inline double G_kalemen(Vec3 L, Vec3 V, Vec3 H, Vec3 N, double alpha) {
	//	return 2.0 * glm::dot(N, L) * glm::dot(N, V) / (1.0 + glm::dot(L, V));
	//}

	inline float BeckmannMicrofacetBRDF_without_F(const glm::vec3 &omega_i, const glm::vec3 &omega_o, const glm::vec3 &omega_h, const glm::vec3 &Ng, float alpha) {
		float d = D_Beckman(Ng, omega_h, alpha);
		float g = G2_height_correlated_beckmann(omega_i, omega_o, omega_h, Ng, alpha);

		float cos_term_wo = glm::dot(Ng, omega_o);
		float cos_term_wi = glm::dot(Ng, omega_i);

		return d * g / (4.0f * cos_term_wo * cos_term_wi);
	}

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

	inline glm::vec3 polar_to_cartesian(float theta, float phi) {
		float sinTheta = sin(theta);
		glm::vec3 v = {
			sinTheta * cos(phi),
			sinTheta * sin(phi),
			cos(theta)
		};
		return v;
	};
}

