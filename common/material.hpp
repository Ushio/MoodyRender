#pragma once

#include <glm/glm.hpp>
#include <glm/ext.hpp>

#include <functional>
#include <strict_variant/variant.hpp>
#include "peseudo_random.hpp"
#include "microfacet.hpp"

namespace rt {
	// z が上, 任意の x, y
	// 一般的な極座標系とも捉えられる
	struct ArbitraryBRDFSpace {
		ArbitraryBRDFSpace(const glm::vec3 &zAxis) {
			zaxis = zAxis;
			if (0.999f < glm::abs(zaxis.z)) {
				xaxis = glm::normalize(glm::cross(glm::vec3(0.0f, 1.0f, 0.0f), zaxis));
			}
			else {
				xaxis = glm::normalize(glm::cross(glm::vec3(0.0f, 0.0f, 1.0f), zaxis));
			}
			yaxis = glm::cross(zaxis, xaxis);
		}
		glm::vec3 localToGlobal(const glm::vec3 v) const {
			return v.x * xaxis + v.y * yaxis + v.z * zaxis;
		}

		// axis on global space
		glm::vec3 xaxis;
		glm::vec3 yaxis;
		glm::vec3 zaxis;
	};

	class LambertianSampler {
	public:
		static glm::vec3 sample(PeseudoRandom *random, const glm::vec3 &Ng) {
			float u1 = random->uniformf();
			float u2 = random->uniformf();
			float r = glm::sqrt(u1);
			float phi = glm::two_pi<float>() * u2;
			glm::vec3 sample(r * glm::cos(phi), r * glm::sin(phi), glm::sqrt(1.0f - u1));
			ArbitraryBRDFSpace space(Ng);
			return space.localToGlobal(sample);
		}
		static float pdf(const glm::vec3 &sampled_wi, const glm::vec3 &Ng) {
			float cosTheta = glm::dot(sampled_wi, Ng);
			if (cosTheta < 0.0f) {
				return 0.0f;
			}
			return cosTheta * glm::one_over_pi<float>();
		}
	};

	class IMaterial {
	public:
		virtual ~IMaterial() {}

		// common member
		glm::vec3 Ng;

		// evaluate emission
		virtual glm::vec3 emission(const glm::vec3 &wo) const {
			return glm::vec3(0.0f);
		}

		// evaluate bxdf
		virtual glm::vec3 bxdf(const glm::vec3 &wo, const glm::vec3 &wi) const = 0;

		// sample wi
		virtual glm::vec3 sample(PeseudoRandom *random, const glm::vec3 &wo) const = 0;

		// pdf for wi
		virtual float pdf(const glm::vec3 &wo, const glm::vec3 &sampled_wi) const = 0;
	};

	class LambertianMaterial : public IMaterial {
	public:
		LambertianMaterial() :Le(0.0f), R(1.0f) {}
		LambertianMaterial(glm::vec3 e, glm::vec3 r) : Le(e), R(r) {}
		glm::vec3 Le;
		glm::vec3 R;

		bool isEmissive() const {
			return glm::any(glm::greaterThanEqual(Le, glm::vec3(glm::epsilon<float>())));
		}
		glm::vec3 emission(const glm::vec3 &wo) const override {
			return Le;
		}
		glm::vec3 bxdf(const glm::vec3 &wo, const glm::vec3 &wi) const override {
			return glm::vec3(R) * glm::one_over_pi<float>();
		}
		glm::vec3 sample(PeseudoRandom *random, const glm::vec3 &wo) const override {
			return LambertianSampler::sample(random, Ng);
		}
		virtual float pdf(const glm::vec3 &wo, const glm::vec3 &sampled_wi) const override {
			return LambertianSampler::pdf(sampled_wi, Ng);
		}
	};

	class SpecularMaterial : public IMaterial {
	public:
		glm::vec3 bxdf(const glm::vec3 &wo, const glm::vec3 &wi) const override {
			return glm::vec3(1.0) * glm::one_over_pi<float>();
		}
		glm::vec3 sample(PeseudoRandom *random, const glm::vec3 &wo) const override {
			return glm::vec3();
		}
		float pdf(const glm::vec3 &wo, const glm::vec3 &sampled_wi) const override {
			return 0.0f;
		}
	};

	class MicrofacetConductorMaterial : public IMaterial {
	public:
		bool useFresnel = true;
		float alpha = 0.2f;

		glm::vec3 bxdf(const glm::vec3 &wo, const glm::vec3 &wi) const override {
			float cos_term_wo = glm::dot(Ng, wo);
			float cos_term_wi = glm::dot(Ng, wi);

			// chi_plus(glm::dot(Ng, omega_i)) * chi_plus(glm::dot(Ng, omega_o))
			if (cos_term_wo <= 0.0f || cos_term_wi <= 0.0f) {
				return glm::vec3();
			}

			glm::vec3 h = glm::normalize(wi + wo);
			float d = D_Beckmann(Ng, h, alpha);
			float g = G2_height_correlated_beckmann(wi, wo, h, Ng, alpha);

			float brdf_without_f = d * g / (4.0f * cos_term_wo * cos_term_wi);

			glm::vec3 brdf = glm::vec3(brdf_without_f);

			if (useFresnel) {
				glm::vec3 eta(0.15557f, 0.42415f, 1.3821f);
				glm::vec3 k(3.6024f, 2.4721f, 1.9155f);

				float cosThetaFresnel = glm::dot(h, wo);
				glm::vec3 f = glm::vec3(
					fresnel_unpolarized(eta.r, k.r, cosThetaFresnel),
					fresnel_unpolarized(eta.g, k.g, cosThetaFresnel),
					fresnel_unpolarized(eta.b, k.b, cosThetaFresnel)
				);
				brdf = f * brdf_without_f;
			}

			return brdf;
		}

		glm::vec3 sample(PeseudoRandom *random, const glm::vec3 &wo) const override {
			return BeckmannImportanceSampler::sample(random, alpha, wo, Ng);
		}
		float pdf(const glm::vec3 &wo, const glm::vec3 &sampled_wi) const override {
			return BeckmannImportanceSampler::pdf(sampled_wi, alpha, wo, Ng);
		}
	};
	class MicrofacetCoupledConductorMaterial : public IMaterial {
	public:
		bool useFresnel = true;
		float alpha = 0.2f;

		glm::vec3 bxdf(const glm::vec3 &wo, const glm::vec3 &wi) const override {
			float cos_term_wo = glm::dot(Ng, wo);
			float cos_term_wi = glm::dot(Ng, wi);

			// chi_plus(glm::dot(Ng, omega_i)) * chi_plus(glm::dot(Ng, omega_o))
			if (cos_term_wo <= 0.0f || cos_term_wi <= 0.0f) {
				return glm::vec3();
			}

			glm::vec3 h = glm::normalize(wi + wo);
			float d = D_Beckmann(Ng, h, alpha);
			float g = G2_height_correlated_beckmann(wi, wo, h, Ng, alpha);

			float brdf_without_f = d * g / (4.0f * cos_term_wo * cos_term_wi);

			glm::vec3 brdf_spec = glm::vec3(brdf_without_f);

			glm::vec3 eta(0.15557f, 0.42415f, 1.3821f);
			glm::vec3 k(3.6024f, 2.4721f, 1.9155f);

			if (useFresnel) {
				float cosThetaFresnel = glm::dot(h, wo);
				glm::vec3 f = glm::vec3(
					fresnel_unpolarized(eta.r, k.r, cosThetaFresnel),
					fresnel_unpolarized(eta.g, k.g, cosThetaFresnel),
					fresnel_unpolarized(eta.b, k.b, cosThetaFresnel)
				);
				brdf_spec = f * brdf_without_f;
			}

			glm::vec3 kLambda = glm::vec3(1.0f);

			if (useFresnel) {
				kLambda = glm::vec3(
					fresnel_unpolarized(eta.r, k.r, 0.0f),
					fresnel_unpolarized(eta.g, k.g, 0.0f),
					fresnel_unpolarized(eta.b, k.b, 0.0f)
				);
			}

			glm::vec3 brdf_diff = kLambda
				* (1.0f - CoupledBRDFConductor::specularAlbedo().sample(alpha, cos_term_wo))
				* (1.0f - CoupledBRDFConductor::specularAlbedo().sample(alpha, cos_term_wi))
				/ (glm::pi<float>() * (1.0f - CoupledBRDFConductor::specularAvgAlbedo().sample(alpha)));

			return brdf_spec + brdf_diff;
		}
		//glm::vec3 sample(PeseudoRandom *random, const glm::vec3 &wo) const override {
		//	// ミックス 重点サンプリング
		//	glm::vec3 wi;
		//	float spAlbedo = CoupledBRDFConductor::specularAlbedo().sample(alpha, glm::dot(Ng, wo));

		//	if (random->uniformf() < spAlbedo) {
		//		wi = BeckmannImportanceSampler::sample(random, alpha, wo, Ng);
		//	}
		//	else {
		//		wi = LambertianSampler::sample(random, Ng);
		//	}
		//	return wi;
		//}
		//float pdf(const glm::vec3 &wo, const glm::vec3 &sampled_wi) const override {
		//	float spAlbedo = CoupledBRDFConductor::specularAlbedo().sample(alpha, glm::dot(Ng, wo));
		//	float pdf_omega = 
		//		spAlbedo * BeckmannImportanceSampler::pdf(sampled_wi, alpha, wo, Ng)
		//		+
		//		(1.0f - spAlbedo) * LambertianSampler::pdf(sampled_wi, Ng);
		//	return pdf_omega;
		//}

		glm::vec3 sample(PeseudoRandom *random, const glm::vec3 &wo) const override {
			glm::vec3 wi;
			float spAlbedo = CoupledBRDFConductor::specularAlbedo().sample(alpha, glm::dot(Ng, wo));

			if (random->uniformf() < spAlbedo) {
				wi = BeckmannImportanceSampler::sample(random, alpha, wo, Ng);
			}
			else {
				float theta = CoupledBRDFConductor::sampler().sampleTheta(alpha, random);
				glm::vec3 sample = polar_to_cartesian(theta, random->uniformf(0.0f, glm::two_pi<float>()));
				ArbitraryBRDFSpace space(Ng);
				return space.localToGlobal(sample);
				// wi = LambertianSampler::sample(random, Ng);
			}
			return wi;
		}
		float pdf(const glm::vec3 &wo, const glm::vec3 &sampled_wi) const override {
			const CoupledBRDFSampler &sampler = CoupledBRDFConductor::sampler();
			float theta = std::acos(glm::dot(Ng, sampled_wi));

			float spAlbedo = CoupledBRDFConductor::specularAlbedo().sample(alpha, glm::dot(Ng, wo));
			float pdf_omega =
				spAlbedo * BeckmannImportanceSampler::pdf(sampled_wi, alpha, wo, Ng)
				+
				(1.0f - spAlbedo) * (1.0f / glm::two_pi<float>()) * (sampler.thetaSize(alpha) * (2.0 / glm::pi<float>())) * sampler.probability(alpha, theta) / std::sin(theta);
				// (1.0f - spAlbedo) * LambertianSampler::pdf(sampled_wi, Ng);
			return pdf_omega;
		}
	};

	typedef strict_variant::variant<LambertianMaterial, SpecularMaterial, MicrofacetConductorMaterial, MicrofacetCoupledConductorMaterial> Material;

	namespace MaterialVisitor {
		struct SetNg {
			SetNg(const glm::vec3 &Ng) : _Ng(Ng) {}
			template <class T>
			void operator()(T &m) {
				m.Ng = _Ng;
			}
			glm::vec3 _Ng;
		};
		struct GetNg {
			template <class T>
			glm::vec3 operator()(T &m) const {
				return m.Ng;
			}
		};
	}
	inline glm::vec3 bxdf_Ng(const Material &m) {
		return strict_variant::apply_visitor(MaterialVisitor::GetNg(), m);
	}
	inline glm::vec3 bxdf_emission(const Material &m, const glm::vec3 &wo) {
		auto f = [](const IMaterial &m, const glm::vec3 &wo) { return m.emission(wo); };
		return strict_variant::apply_visitor(std::bind(f, std::placeholders::_1, wo), m);
	}
	inline glm::vec3 bxdf_evaluate(const Material &m, const glm::vec3 &wo, const glm::vec3 &wi) {
		auto f = [](const IMaterial &m, const glm::vec3 &wo, const glm::vec3 &wi) { return m.bxdf(wo, wi); };
		return strict_variant::apply_visitor(std::bind(f, std::placeholders::_1, wo, wi), m);
	}
	inline glm::vec3 bxdf_sample(const Material &m, PeseudoRandom *random, const glm::vec3 &wo) {
		auto f = [](const IMaterial &m, PeseudoRandom *random, const glm::vec3 &wo) { return m.sample(random, wo); };
		return strict_variant::apply_visitor(std::bind(f, std::placeholders::_1, random, wo), m);
	}
	inline float bxdf_pdf(const Material &m, const glm::vec3 &wo, const glm::vec3 &sampled_wi) {
		auto f = [](const IMaterial &m, const glm::vec3 &wo, const glm::vec3 &sampled_wi) { return m.pdf(wo, sampled_wi); };
		return strict_variant::apply_visitor(std::bind(f, std::placeholders::_1, wo, sampled_wi), m);
	}
}