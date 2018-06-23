#pragma once

#include <glm/glm.hpp>
#include <glm/ext.hpp>

#include <functional>
#include <strict_variant/variant.hpp>
#include "peseudo_random.hpp"
#include "microfacet.hpp"
#include "coordinate.hpp"

namespace rt {
	class LambertianSampler {
	public:
		static glm::dvec3 sample(PeseudoRandom *random, const glm::dvec3 &Ng) {
			double u1 = random->uniform();
			double u2 = random->uniform();
			double r = glm::sqrt(u1);
			double phi = glm::two_pi<double>() * u2;
			glm::dvec3 sample(r * glm::cos(phi), r * glm::sin(phi), glm::sqrt(1.0 - u1));
			ArbitraryBRDFSpace space(Ng);
			return space.localToGlobal(sample);
		}
		static double pdf(const glm::dvec3 &sampled_wi, const glm::dvec3 &Ng) {
			double cosTheta = glm::dot(sampled_wi, Ng);
			if (cosTheta < 0.0) {
				return 0.0;
			}
			return cosTheta * glm::one_over_pi<double>();
		}
	};

	class IMaterial {
	public:
		virtual ~IMaterial() {}

		// common member
		glm::dvec3 Ng;

		virtual bool isEmission() const {
			return false;
		}

		// evaluate emission
		virtual glm::dvec3 emission(const glm::dvec3 &wo) const {
			return glm::dvec3(0.0);
		}

		// evaluate bxdf
		virtual glm::dvec3 bxdf(const glm::dvec3 &wo, const glm::dvec3 &wi) const = 0;

		// sample wi
		virtual glm::dvec3 sample(PeseudoRandom *random, const glm::dvec3 &wo) const = 0;

		// pdf for wi
		virtual double pdf(const glm::dvec3 &wo, const glm::dvec3 &sampled_wi) const = 0;
	};

	class LambertianMaterial : public IMaterial {
	public:
		LambertianMaterial() :Le(0.0), R(1.0) {}
		LambertianMaterial(glm::dvec3 e, glm::dvec3 r) : Le(e), R(r) {}
		glm::dvec3 Le;
		glm::dvec3 R;

		bool isEmission() const override {
			return glm::any(glm::greaterThanEqual(Le, glm::dvec3(glm::epsilon<double>())));
		}
		glm::dvec3 emission(const glm::dvec3 &wo) const override {
			return Le;
		}
		glm::dvec3 bxdf(const glm::dvec3 &wo, const glm::dvec3 &wi) const override {
			if (glm::dot(Ng, wi) < 0.0 || glm::dot(Ng, wo) < 0.0) {
				return glm::dvec3(0.0);
			}
			return glm::dvec3(R) * glm::one_over_pi<double>();
		}
		glm::dvec3 sample(PeseudoRandom *random, const glm::dvec3 &wo) const override {
			return LambertianSampler::sample(random, Ng);
		}
		virtual double pdf(const glm::dvec3 &wo, const glm::dvec3 &sampled_wi) const override {
			return LambertianSampler::pdf(sampled_wi, Ng);
		}
	};

	class SpecularMaterial : public IMaterial {
	public:
		glm::dvec3 bxdf(const glm::dvec3 &wo, const glm::dvec3 &wi) const override {
			return glm::dvec3(1.0) * glm::one_over_pi<double>();
		}
		glm::dvec3 sample(PeseudoRandom *random, const glm::dvec3 &wo) const override {
			return glm::dvec3();
		}
		double pdf(const glm::dvec3 &wo, const glm::dvec3 &sampled_wi) const override {
			return 0.0;
		}
	};

	class MicrofacetConductorMaterial : public IMaterial {
	public:
		bool useFresnel = true;
		double alpha = 0.2;

		glm::dvec3 bxdf(const glm::dvec3 &wo, const glm::dvec3 &wi) const override {
			double cos_term_wo = glm::dot(Ng, wo);
			double cos_term_wi = glm::dot(Ng, wi);

			// chi_plus(glm::dot(Ng, omega_i)) * chi_plus(glm::dot(Ng, omega_o))
			if (cos_term_wo <= 0.0 || cos_term_wi <= 0.0) {
				return glm::dvec3();
			}

			glm::dvec3 h = glm::normalize(wi + wo);
			double d = D_Beckmann(Ng, h, alpha);
			// double g = G2_height_correlated_beckmann(wi, wo, h, Ng, alpha);
			double g = G2_v_cavity(wi, wo, h, Ng);

			double brdf_without_f = d * g / (4.0 * cos_term_wo * cos_term_wi);

			glm::dvec3 brdf = glm::dvec3(brdf_without_f);

			if (useFresnel) {
				glm::dvec3 eta(0.15557, 0.42415, 1.3821);
				glm::dvec3 k(3.6024, 2.4721, 1.9155);

				double cosThetaFresnel = glm::dot(h, wo);
				glm::dvec3 f = glm::dvec3(
					fresnel_unpolarized(eta.r, k.r, cosThetaFresnel),
					fresnel_unpolarized(eta.g, k.g, cosThetaFresnel),
					fresnel_unpolarized(eta.b, k.b, cosThetaFresnel)
				);
				brdf = f * brdf_without_f;
			}

			return brdf;
		}
		glm::dvec3 sample(PeseudoRandom *random, const glm::dvec3 &wo) const override {
			return VCavityBeckmannVisibleNormalSampler::sample(random, alpha, wo, Ng);
		}
		double pdf(const glm::dvec3 &wo, const glm::dvec3 &sampled_wi) const override {
			return VCavityBeckmannVisibleNormalSampler::pdf(sampled_wi, alpha, wo, Ng);
		}
		
		//glm::dvec3 sample(PeseudoRandom *random, const glm::dvec3 &wo) const override {
		//	return BeckmannImportanceSampler::sample(random, alpha, wo, Ng);
		//}
		//double pdf(const glm::dvec3 &wo, const glm::dvec3 &sampled_wi) const override {
		//	return BeckmannImportanceSampler::pdf(sampled_wi, alpha, wo, Ng);
		//}
	};
	class MicrofacetCoupledConductorMaterial : public IMaterial {
	public:
		bool useFresnel = true;
		double alpha = 0.2;

		glm::dvec3 bxdf(const glm::dvec3 &wo, const glm::dvec3 &wi) const override {
			double cos_term_wo = glm::dot(Ng, wo);
			double cos_term_wi = glm::dot(Ng, wi);

			// chi_plus(glm::dot(Ng, omega_i)) * chi_plus(glm::dot(Ng, omega_o))
			if (cos_term_wo <= 0.0 || cos_term_wi <= 0.0) {
				return glm::dvec3();
			}

			glm::dvec3 h = glm::normalize(wi + wo);
			double d = D_Beckmann(Ng, h, alpha);
			// double g = G2_height_correlated_beckmann(wi, wo, h, Ng, alpha);
			double g = G2_v_cavity(wi, wo, h, Ng);

			double brdf_without_f = d * g / (4.0 * cos_term_wo * cos_term_wi);

			glm::dvec3 brdf_spec = glm::dvec3(brdf_without_f);

			glm::dvec3 eta(0.15557, 0.42415, 1.3821);
			glm::dvec3 k(3.6024, 2.4721, 1.9155);

			if (useFresnel) {
				double cosThetaFresnel = glm::dot(h, wo);
				glm::dvec3 f = glm::dvec3(
					fresnel_unpolarized(eta.r, k.r, cosThetaFresnel),
					fresnel_unpolarized(eta.g, k.g, cosThetaFresnel),
					fresnel_unpolarized(eta.b, k.b, cosThetaFresnel)
				);
				brdf_spec = f * brdf_without_f;
			}

			glm::dvec3 kLambda = glm::dvec3(1.0);

			if (useFresnel) {
				kLambda = glm::dvec3(
					fresnel_unpolarized(eta.r, k.r, 0.0),
					fresnel_unpolarized(eta.g, k.g, 0.0),
					fresnel_unpolarized(eta.b, k.b, 0.0)
				);
			}

			glm::dvec3 brdf_diff = kLambda
				* (1.0 - CoupledBRDFConductor::specularAlbedo().sample(alpha, cos_term_wo))
				* (1.0 - CoupledBRDFConductor::specularAlbedo().sample(alpha, cos_term_wi))
				/ (glm::pi<double>() * (1.0 - CoupledBRDFConductor::specularAvgAlbedo().sample(alpha)));

			return brdf_spec + brdf_diff;
		}

		glm::dvec3 sample(PeseudoRandom *random, const glm::dvec3 &wo) const override {
			glm::dvec3 wi;
			double spAlbedo = CoupledBRDFConductor::specularAlbedo().sample(alpha, glm::dot(Ng, wo));

			if (random->uniform() < spAlbedo) {
				wi = VCavityBeckmannVisibleNormalSampler::sample(random, alpha, wo, Ng);
			}
			else {
				double theta = CoupledBRDFConductor::sampler().sampleTheta(alpha, random);
				glm::dvec3 sample = polar_to_cartesian(theta, random->uniform(0.0, glm::two_pi<double>()));
				ArbitraryBRDFSpace space(Ng);
				return space.localToGlobal(sample);
				// wi = LambertianSampler::sample(random, Ng);
			}
			return wi;
		}
		double pdf(const glm::dvec3 &wo, const glm::dvec3 &sampled_wi) const override {
			const CoupledBRDFSampler &sampler = CoupledBRDFConductor::sampler();
			double theta = std::acos(glm::dot(Ng, sampled_wi));

			double spAlbedo = CoupledBRDFConductor::specularAlbedo().sample(alpha, glm::dot(Ng, wo));
			double pdf_omega =
				spAlbedo * VCavityBeckmannVisibleNormalSampler::pdf(sampled_wi, alpha, wo, Ng)
				+
				(1.0 - spAlbedo) * (1.0 / glm::two_pi<double>()) * (sampler.thetaSize(alpha) * (2.0 / glm::pi<double>())) * sampler.probability(alpha, theta) / std::sin(theta);
				// (1.0 - spAlbedo) * LambertianSampler::pdf(sampled_wi, Ng);
			return pdf_omega;
		}
	};
	
	class MicrofacetCoupledDielectricsMaterial : public IMaterial {
	public:
		double alpha = 0.2;
		glm::dvec3 Cd = glm::dvec3(1.0);

		glm::dvec3 bxdf(const glm::dvec3 &wo, const glm::dvec3 &wi) const override {
			double cos_term_wo = glm::dot(Ng, wo);
			double cos_term_wi = glm::dot(Ng, wi);

			// chi_plus(glm::dot(Ng, omega_i)) * chi_plus(glm::dot(Ng, omega_o))
			if (cos_term_wo <= 0.0 || cos_term_wi <= 0.0) {
				return glm::dvec3();
			}

			glm::dvec3 h = glm::normalize(wi + wo);
			double d = D_Beckmann(Ng, h, alpha);
			// double g = G2_height_correlated_beckmann(wi, wo, h, Ng, alpha);
			double g = G2_v_cavity(wi, wo, h, Ng);

			double brdf_without_f = d * g / (4.0 * cos_term_wo * cos_term_wi);

			glm::dvec3 brdf_spec = glm::dvec3(brdf_without_f);

			{
				double cosThetaFresnel = glm::dot(h, wo);
				glm::dvec3 f = glm::dvec3(fresnel_dielectrics(cosThetaFresnel));
				brdf_spec = f * brdf_without_f;
			}

			glm::dvec3 kLambda = Cd;

			glm::dvec3 brdf_diff = kLambda
				* (1.0 - CoupledBRDFDielectrics::specularAlbedo().sample(alpha, cos_term_wo))
				* (1.0 - CoupledBRDFDielectrics::specularAlbedo().sample(alpha, cos_term_wi))
				/ (glm::pi<double>() * (1.0 - CoupledBRDFDielectrics::specularAvgAlbedo().sample(alpha)));

			return brdf_spec + brdf_diff;
		}

		glm::dvec3 sample(PeseudoRandom *random, const glm::dvec3 &wo) const override {
			glm::dvec3 wi;
			double spAlbedo = CoupledBRDFDielectrics::specularAlbedo().sample(alpha, glm::dot(Ng, wo));

			if (random->uniform() < spAlbedo) {
				wi = VCavityBeckmannVisibleNormalSampler::sample(random, alpha, wo, Ng);
			}
			else {
				double theta = CoupledBRDFDielectrics::sampler().sampleTheta(alpha, random);
				glm::dvec3 sample = polar_to_cartesian(theta, random->uniform(0.0, glm::two_pi<double>()));
				ArbitraryBRDFSpace space(Ng);
				return space.localToGlobal(sample);
				// wi = LambertianSampler::sample(random, Ng);
			}
			return wi;
		}
		double pdf(const glm::dvec3 &wo, const glm::dvec3 &sampled_wi) const override {
			const CoupledBRDFSampler &sampler = CoupledBRDFDielectrics::sampler();
			double theta = std::acos(glm::dot(Ng, sampled_wi));

			double spAlbedo = CoupledBRDFDielectrics::specularAlbedo().sample(alpha, glm::dot(Ng, wo));
			double pdf_omega =
				spAlbedo * VCavityBeckmannVisibleNormalSampler::pdf(sampled_wi, alpha, wo, Ng)
				+
				(1.0 - spAlbedo) * (1.0 / glm::two_pi<double>()) * (sampler.thetaSize(alpha) * (2.0 / glm::pi<double>())) * sampler.probability(alpha, theta) / std::sin(theta);
			// (1.0 - spAlbedo) * LambertianSampler::pdf(sampled_wi, Ng);
			return pdf_omega;
		}
	};

	class Material {
	public:
		typedef strict_variant::variant<
			LambertianMaterial,
			SpecularMaterial,
			MicrofacetConductorMaterial,
			MicrofacetCoupledConductorMaterial,
			MicrofacetCoupledDielectricsMaterial
		> MaterialType;

		Material() {}
		Material(const MaterialType &m) :_material(m) { }
		Material &operator=(const MaterialType &m) { _material = m; return *this; }

		IMaterial *operator->() {
			return strict_variant::apply_visitor([](IMaterial &m) { return &m; }, _material);
		}
		const IMaterial *operator->() const {
			return strict_variant::apply_visitor([](const IMaterial &m) { return &m; }, _material);
		}
		
		template <class T>
		T *get() {
			return _material.get<T>();
		}
		template <class T>
		const T *get() const {
			return _material.get<T>();
		}
	private:
		MaterialType _material;
	};
}