#pragma once

#include <glm/glm.hpp>
#include <glm/ext.hpp>

namespace rt {

#define MICROCYLINDER_ASSERT_ENABLE 1

#if MICROCYLINDER_ASSERT_ENABLE
#define MICROCYLINDER_ASSERT(expect_true) if((expect_true) == 0) { __debugbreak(); }
#else
#define MICROCYLINDER_ASSERT(expect_true) ;
#endif

	inline double normalized_gaussian(double beta, double theta) {
		return std::exp(-theta * theta / (2.0 * beta * beta)) / (std::sqrt(glm::two_pi<double>() * beta * beta));
	}

	// ex)
	// air to glass
	// eta_t(glass) = 1.5, eta_i(air) = 1.0
	inline double microcylinder_fresnel_dielectrics(double cosTheta, double eta_t, double eta_i) {
		auto sqr = [](double x) { return x * x; };

		double c = cosTheta;
		double g = std::sqrt(sqr(eta_t) / sqr(eta_i) - 1.0 + sqr(c));

		double a = 0.5 * sqr(g - c) / sqr(g + c);
		double b = 1.0 + sqr(c * (g + c) - 1.0) / sqr(c * (g - c) + 1.0);
		return a * b;
	}

	inline double fr_cosTheta(double theta_d, double phi_d) {
		return std::cos(theta_d) * std::cos(phi_d * 0.5);
	}

	inline double scattering_rs(double phi_d, double theta_h, double gamma_s) {
		return std::cos(phi_d * 0.5) * normalized_gaussian(gamma_s, theta_h);
	}
	inline double scattering_rv(double theta_h, double gamma_v, double kd, double cosThetaI, double cosThetaO) {
		return ((1.0 - kd) * normalized_gaussian(gamma_v, theta_h) + kd) / (cosThetaI + cosThetaO);
	}

	struct ThreadParameters {
		double gamma_s;
		double gamma_v;
		double kd;
		double eta_t;
		glm::dvec3 A;
	};

	inline glm::dvec3 microcylinder_bsdf(double theta_d, double theta_h, double phi_d, double cosThetaI, double cosThetaO, ThreadParameters params) {
		double eta_i = 1.0;
		MICROCYLINDER_ASSERT(1.0 <= eta_i);

		// f_r,s
		// return glm::dvec3(1.0) * scattering_rs(phi_d, theta_h, gamma_s) / std::pow(cos(theta_d), 2);

		double Fr_cosTheta_i = fr_cosTheta(theta_d, phi_d);
		auto safeSqrt = [](double x) {
			return std::sqrt(std::max(x, 0.0));
		};

		double Fr = microcylinder_fresnel_dielectrics(Fr_cosTheta_i, params.eta_t, eta_i);
		double Ft = (1.0 - Fr);
		double F = Ft * Ft;

		double rs = Fr * scattering_rs(phi_d, theta_h, params.gamma_s);
		double rv = F * scattering_rv(theta_h, params.gamma_v, params.kd, cosThetaI, cosThetaO);

		// f_r,v
		// return (scattering_rv(theta_h, gamma_v, kd, *cosThetaI, cosThetaO) * A) / std::pow(cos(theta_d), 2);

		return (glm::dvec3(rs) + rv * glm::dvec3(params.A)) / std::pow(cos(theta_d), 2);
	}

	struct ThreadGeometricParameters {
		double theta_d;
		double theta_h;
		double phi_d;
		double cosPhiI;
		double cosPhiO;
		double cosThetaI;
		double cosThetaO;
		double psi_d;
		double cosPsiI;
		double cosPsiO;
	};

	inline bool parameterize(glm::dvec3 u, glm::dvec3 v, glm::dvec3 n, glm::dvec3 wi, glm::dvec3 wo, ThreadGeometricParameters *params) {
		MICROCYLINDER_ASSERT(std::abs(glm::length2(u) - 1.0) < 1.0e-4);
		MICROCYLINDER_ASSERT(std::abs(glm::length2(n) - 1.0) < 1.0e-4);
		MICROCYLINDER_ASSERT(std::abs(glm::length2(wi) - 1.0) < 1.0e-4);
		MICROCYLINDER_ASSERT(std::abs(glm::length2(wo) - 1.0) < 1.0e-4);

		double sinThetaI = glm::clamp(glm::dot(wi, u), -1.0, 1.0); // clamp for asin stability
		double thetaI = std::asin(sinThetaI);
		double sinThetaO = glm::clamp(glm::dot(wo, u), -1.0, 1.0); // clamp for asin stability
		double thetaO = std::asin(sinThetaO);

		glm::dvec3 wi_on_normal = glm::normalize(wi - u * sinThetaI);
		glm::dvec3 wo_on_normal = glm::normalize(wo - u * sinThetaO);

		// Φが定義できない
		if (glm::all(glm::isfinite(wi_on_normal)) == false) {
			return false;
		}
		if (glm::all(glm::isfinite(wo_on_normal)) == false) {
			return false;
		}

		double cosPhiD = glm::clamp(glm::dot(wi_on_normal, wo_on_normal), -1.0, 1.0); // clamp for acos stability

		double phi_d = std::acos(cosPhiD);
		double theta_h = (thetaI + thetaO) * 0.5;
		double theta_d = (thetaI - thetaO) * 0.5;

		double cosThetaO = std::cos(thetaO);

		glm::dvec3 wi_on_tangent_normal = glm::normalize(wi - v * glm::dot(wi, v));
		glm::dvec3 wo_on_tangent_normal = glm::normalize(wo - v * glm::dot(wo, v));

		// ψが定義できない
		if (glm::all(glm::isfinite(wi_on_tangent_normal)) == false) {
			return false;
		}
		if (glm::all(glm::isfinite(wo_on_tangent_normal)) == false) {
			return false;
		}
		double cosPsiD = glm::clamp(glm::dot(wi_on_tangent_normal, wo_on_tangent_normal), -1.0, 1.0); // clamp for acos stability
		double psi_d = std::acos(cosPsiD);

		params->theta_d = theta_d;
		params->theta_h = theta_h;
		params->cosPhiI = glm::dot(n, wi_on_normal);
		params->cosPhiO = glm::dot(n, wo_on_normal);
		params->phi_d = phi_d;
		params->cosThetaI = std::cos(thetaI);
		params->cosThetaO = std::cos(thetaO);
		params->psi_d = psi_d;
		params->cosPsiI = glm::dot(n, wi_on_tangent_normal);
		params->cosPsiO = glm::dot(n, wo_on_tangent_normal);

		return true;
	}

	inline glm::dvec3 microcylinder_bsdf(glm::dvec3 u, glm::dvec3 v, glm::dvec3 n, glm::dvec3 wi, glm::dvec3 wo, ThreadParameters params, double *cosThetaI) {
		MICROCYLINDER_ASSERT(0.0 <= params.gamma_s);
		MICROCYLINDER_ASSERT(0.0 <= params.gamma_v);
		MICROCYLINDER_ASSERT(1.0 <= params.eta_t);
		MICROCYLINDER_ASSERT(0.0 <= params.kd && params.kd <= 1.0);
		MICROCYLINDER_ASSERT(glm::all(glm::greaterThanEqual(params.A, glm::dvec3(0.0))));

		ThreadGeometricParameters p;
		if (parameterize(u, v, n, wi, wo, &p) == false) {
			return glm::dvec3(0.0);
		}

		*cosThetaI = p.cosThetaI;
		MICROCYLINDER_ASSERT(0.0 <= *cosThetaI);

		return microcylinder_bsdf(p.theta_d, p.theta_h, p.phi_d, p.cosThetaI, p.cosThetaO, params);
	}

	inline double u_gaussian(double x) {
		static const double sd = glm::radians(20.0);
		static const double sqr_sd_2 = sd * sd * 2.0;
		return std::exp(-x * x / sqr_sd_2);
	}

	inline double microcylinder_M(double cosPhiI, double cosPhiO, double phi_d) {
		double m_i = std::max(cosPhiI, 0.0);
		double m_o = std::max(cosPhiO, 0.0);
		double u = u_gaussian(phi_d);
		double corrated = std::min(m_i, m_o);
		double uncorrated = m_i * m_o;
		return glm::mix(uncorrated, corrated, u);
	}
	inline double microcylinder_P(double cosPsiI, double cosPsiO, double psi_d) {
		double m_i = std::max(cosPsiI, 0.0);
		double m_o = std::max(cosPsiO, 0.0);
		double u = u_gaussian(psi_d);
		double corrated = std::min(m_i, m_o);
		double uncorrated = m_i * m_o;
		return glm::mix(uncorrated, corrated, u);
	}

	inline glm::dvec3 brdf_x_cosTheta(glm::dvec3 u, glm::dvec3 v, glm::dvec3 n, glm::dvec3 wi, glm::dvec3 wo, double alpha_0, double alpha_1, const std::vector<double> &tangent_offsets_u, const std::vector<double> &tangent_offsets_v, ThreadParameters params_0, ThreadParameters params_1, double *cosThetaI) {
		MICROCYLINDER_ASSERT(0.0 <= alpha_0 && alpha_0 <= 1.0);
		MICROCYLINDER_ASSERT(0.0 <= alpha_1 && alpha_1 <= 1.0);

		MICROCYLINDER_ASSERT(glm::abs(glm::dot(u, v)) < 1.0e-3);
		MICROCYLINDER_ASSERT(glm::abs(glm::dot(v, n)) < 1.0e-3);
		MICROCYLINDER_ASSERT(glm::abs(glm::dot(n, v)) < 1.0e-3);

		int Nu = tangent_offsets_u.size();
		int Nv = tangent_offsets_v.size();

		glm::dvec3 uValue;
		glm::dvec3 vValue;

		double Q = 0.0;

		for (double tangent_offset : tangent_offsets_u) {
			glm::dvec3 v_shade = v;
			glm::dvec3 u_shade = glm::rotate(u, glm::radians(tangent_offset), v_shade);
			glm::dvec3 n_shade = glm::rotate(n, glm::radians(tangent_offset), v_shade);
			MICROCYLINDER_ASSERT(std::abs(glm::dot(u_shade, u) - std::cos(glm::radians(tangent_offset))) < 0.0001);
			MICROCYLINDER_ASSERT(std::abs(glm::dot(n_shade, n) - std::cos(glm::radians(tangent_offset))) < 0.0001);

			ThreadGeometricParameters p;
			if (parameterize(u_shade, v_shade, n_shade, wi, wo, &p) == false) {
				continue;
			}
			glm::dvec3 fr = microcylinder_bsdf(p.theta_d, p.theta_h, p.phi_d, p.cosThetaI, p.cosThetaO, params_0);

			// simple
			// uValue += fr * p.cosThetaI;

			// shadowing masking
			// double m = microcylinder_M(p.cosPhiI, p.cosPhiO, p.phi_d);
			// uValue += m * fr * p.cosThetaI;

			// shadowing masking + reweighting
			double m_value = microcylinder_M(p.cosPhiI, p.cosPhiO, p.phi_d);
			double p_value = microcylinder_P(p.cosPsiI, p.cosPsiO, p.psi_d);
			uValue += p_value * m_value * fr * p.cosThetaI;
			Q += alpha_0 * p_value / Nu;
		}
		uValue /= Nu;

		for (double tangent_offset : tangent_offsets_v) {
			glm::dvec3 v_shade = u;
			glm::dvec3 u_shade = glm::rotate(v, glm::radians(tangent_offset), v_shade);
			glm::dvec3 n_shade = glm::rotate(n, glm::radians(tangent_offset), v_shade);

			MICROCYLINDER_ASSERT(std::abs(glm::dot(u_shade, v) - std::cos(glm::radians(tangent_offset))) < 0.0001);
			MICROCYLINDER_ASSERT(std::abs(glm::dot(n_shade, n) - std::cos(glm::radians(tangent_offset))) < 0.0001);

			ThreadGeometricParameters p;
			if (parameterize(u_shade, v_shade, n_shade, wi, wo, &p) == false) {
				continue;
			}
			glm::dvec3 fr = microcylinder_bsdf(p.theta_d, p.theta_h, p.phi_d, p.cosThetaI, p.cosThetaO, params_1);

			// simple
			// vValue += fr * p.cosThetaI;

			// double m = microcylinder_M(p.cosPhiI, p.cosPhiO, p.phi_d);
			// vValue += m * fr * p.cosThetaI;

			// shadowing masking + reweighting
			double m_value = microcylinder_M(p.cosPhiI, p.cosPhiO, p.phi_d);
			double p_value = microcylinder_P(p.cosPsiI, p.cosPsiO, p.psi_d);
			vValue += p_value * m_value * fr * p.cosThetaI;
			Q += alpha_1 * p_value / Nv;
		}
		vValue /= Nv;

		*cosThetaI = glm::abs(glm::dot(n, wi));
		glm::dvec3 fr_cosTheta = uValue * alpha_0 + vValue * alpha_1;

		if (0.0 < Q) {
			Q += (1.0 - alpha_0 - alpha_1) * glm::dot(wi, n);

			fr_cosTheta /= Q;
		}

		return fr_cosTheta;
	}

	inline ThreadParameters a_both() {
		ThreadParameters p;
		p.gamma_s = glm::radians(12.0);
		p.gamma_v = glm::radians(24.0);
		p.eta_t = 1.46;
		p.kd = 0.3;
		p.A = glm::vec3(0.2, 0.8, 1.0) * 0.3;
		return p;
	}

	inline ThreadParameters b_flat() {
		ThreadParameters p;
		p.gamma_s = glm::radians(5.0);
		p.gamma_v = glm::radians(10.0);
		p.eta_t = 1.345;
		p.kd = 0.2;
		p.A = glm::vec3(1.0, 0.95, 0.05) * 0.12;
		return p;
	}

	inline ThreadParameters b_twisted() {
		ThreadParameters p;
		p.gamma_s = glm::radians(18.0);
		p.gamma_v = glm::radians(32.0);
		p.eta_t = 1.345;
		p.kd = 0.3;
		p.A = glm::vec3(1.0, 0.95, 0.05) * 0.16;
		return p;
	}

	inline ThreadParameters c_flat() {
		ThreadParameters p;
		p.gamma_s = glm::radians(2.5);
		p.gamma_v = glm::radians(5.0);
		p.eta_t = 1.539;
		p.kd = 0.1;
		p.A = glm::vec3(1.0, 0.37, 0.3) * 0.035;
		return p;
	}

	inline ThreadParameters c_twisted() {
		ThreadParameters p;
		p.gamma_s = glm::radians(30.0);
		p.gamma_v = glm::radians(60.0);
		p.eta_t = 1.539;
		p.kd = 0.7;
		p.A = glm::vec3(1.0, 0.37, 0.3) * 0.2;
		return p;
	}

	inline ThreadParameters e_dir1() {
		ThreadParameters p;
		p.gamma_s = glm::radians(4.0);
		p.gamma_v = glm::radians(8.0);
		p.eta_t = 1.345;
		p.kd = 0.1;
		p.A = glm::vec3(0.1, 1.0, 0.4) * 0.2;
		return p;
	}
	inline ThreadParameters e_dir2() {
		ThreadParameters p;
		p.gamma_s = glm::radians(5.0);
		p.gamma_v = glm::radians(10.0);
		p.eta_t = 1.345;
		p.kd = 0.1;
		p.A = glm::vec3(1.0, 0.0, 0.1) * 0.6;
		return p;
	}

	inline ThreadParameters f_velvet() {
		ThreadParameters p;
		p.gamma_s = glm::radians(6.0);
		p.gamma_v = glm::radians(12.0);
		p.eta_t = 1.46;
		p.kd = 0.1;
		p.A = glm::vec3(0.05, 0.02, 0.0) * 0.3;
		return p;
	}

	inline void microcylinder_validation() {
		MICROCYLINDER_ASSERT(std::abs(u_gaussian(0.0) - 1.0) < 1.0e-5);
		MICROCYLINDER_ASSERT(std::abs(microcylinder_M(0.5, 0.5, 0.0) - 0.5) < 1.0e-5);
		MICROCYLINDER_ASSERT(std::abs(microcylinder_P(0.5, 0.5, 0.0) - 0.5) < 1.0e-5);
		MICROCYLINDER_ASSERT(std::abs(microcylinder_M(0.5, 0.5, 10000.0) - 0.25) < 1.0e-5);
		MICROCYLINDER_ASSERT(std::abs(microcylinder_P(0.5, 0.5, 10000.0) - 0.25) < 1.0e-5);

		{
			rt::XoroshiroPlus128 random;
			for (int i = 0; i < 100000; ++i) {
				double x = random.uniform(0.0, 10.0);
				double v1 = u_gaussian( x);
				double v2 = u_gaussian(-x);
				MICROCYLINDER_ASSERT(0.0 <= v1 && v1 <= 1.0);
				MICROCYLINDER_ASSERT(0.0 <= v2 && v2 <= 1.0);
				MICROCYLINDER_ASSERT(std::abs(v1 - v2) < 1.0e-5);
			}
		}
		
		glm::dvec3 u = glm::dvec3(0, 0, 1);
		glm::dvec3 v = glm::dvec3(1, 0, 0);
		glm::dvec3 n = glm::dvec3(0, 1, 0);

		// [cosThetaI, cosThetaO, theta_h, theta_d]
		{
			// theta = -45deg, 45deg
			glm::dvec3 wi = glm::normalize(glm::dvec3(0, 1, +1));
			glm::dvec3 wo = glm::normalize(glm::dvec3(0, 1, -1));

			rt::XoroshiroPlus128 random;

			for (int i = 0; i < 10000; ++i) {
				glm::dvec3 any_wi = glm::rotate(wi, random.uniform(0.0, glm::two_pi<double>()), u);
				glm::dvec3 any_wo = glm::rotate(wo, random.uniform(0.0, glm::two_pi<double>()), u);
				ThreadGeometricParameters p;
				if (parameterize(u, v, n, any_wi, any_wo, &p) == false) {
					MICROCYLINDER_ASSERT(0);
				}

				MICROCYLINDER_ASSERT(std::abs(p.cosThetaI - glm::cos(glm::radians(45.0))) < 1.0e-5);
				MICROCYLINDER_ASSERT(std::abs(p.cosThetaO - glm::cos(glm::radians(45.0))) < 1.0e-5);
				MICROCYLINDER_ASSERT(std::abs(p.theta_h) < 1.0e-5);
				MICROCYLINDER_ASSERT(std::abs(p.theta_d - glm::radians(45.0)) < 1.0e-5);
			}
		}

		// [cosPhiI, cosPhiO, phi_d]
		{
			// phi = -45deg, 45deg
			glm::dvec3 wi = glm::normalize(glm::dvec3(+1, 1, 0));
			glm::dvec3 wo = glm::normalize(glm::dvec3(-1, 1, 0));

			rt::XoroshiroPlus128 random;

			for (int i = 0; i < 10000; ++i) {
				glm::dvec3 axis = glm::normalize(glm::cross(wi, u));
				glm::dvec3 any_wi = glm::rotate(wi, random.uniform(-glm::pi<double>() * 0.5, glm::pi<double>() * 0.5), axis);
				ThreadGeometricParameters p;
				if (parameterize(u, v, n, any_wi, wo, &p) == false) {
					MICROCYLINDER_ASSERT(0);
				}
				MICROCYLINDER_ASSERT(std::abs(p.cosPhiI - glm::cos(glm::radians(45.0))) < 1.0e-5);
				MICROCYLINDER_ASSERT(std::abs(p.phi_d - glm::radians(90.0)) < 1.0e-5);
			}
			for (int i = 0; i < 10000; ++i) {
				glm::dvec3 axis = glm::normalize(glm::cross(wo, u));
				glm::dvec3 any_wo = glm::rotate(wo, random.uniform(-glm::pi<double>() * 0.5, glm::pi<double>() * 0.5), axis);
				ThreadGeometricParameters p;
				if (parameterize(u, v, n, wi, any_wo, &p) == false) {
					MICROCYLINDER_ASSERT(0);
				}
				MICROCYLINDER_ASSERT(std::abs(p.cosPhiO - glm::cos(glm::radians(45.0))) < 1.0e-5);
				MICROCYLINDER_ASSERT(std::abs(p.phi_d - glm::radians(90.0)) < 1.0e-5);
			}
		}

		// [cosPsiI, cosPsiO, psi_d]
		{
			// psi = -45deg, 45deg
			glm::dvec3 wi = glm::normalize(glm::dvec3(0, 1, +1));
			glm::dvec3 wo = glm::normalize(glm::dvec3(0, 1, -1));

			rt::XoroshiroPlus128 random;

			for (int i = 0; i < 10000; ++i) {
				glm::dvec3 axis = glm::normalize(glm::cross(wi, v));
				glm::dvec3 any_wi = glm::rotate(wi, random.uniform(-glm::pi<double>() * 0.5, glm::pi<double>() * 0.5), axis);
				ThreadGeometricParameters p;
				if (parameterize(u, v, n, any_wi, wo, &p) == false) {
					MICROCYLINDER_ASSERT(0);
				}
				MICROCYLINDER_ASSERT(std::abs(p.cosPsiI - glm::cos(glm::radians(45.0))) < 1.0e-5);
				MICROCYLINDER_ASSERT(std::abs(p.psi_d - glm::radians(90.0)) < 1.0e-5);
			}
			for (int i = 0; i < 10000; ++i) {
				glm::dvec3 axis = glm::normalize(glm::cross(wo, v));
				glm::dvec3 any_wo = glm::rotate(wo, random.uniform(-glm::pi<double>() * 0.5, glm::pi<double>() * 0.5), axis);
				ThreadGeometricParameters p;
				if (parameterize(u, v, n, wi, any_wo, &p) == false) {
					MICROCYLINDER_ASSERT(0);
				}
				MICROCYLINDER_ASSERT(std::abs(p.cosPsiO - glm::cos(glm::radians(45.0))) < 1.0e-5);
				MICROCYLINDER_ASSERT(std::abs(p.psi_d - glm::radians(90.0)) < 1.0e-5);
			}
		}
	}
}