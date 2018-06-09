#include "ofApp.h"
#include "serializable_buffer.hpp"
#include "peseudo_random.hpp"
#include "microfacet.hpp"
#include "online.hpp"
#include "stopwatch.hpp"
#include "composite_simpson.hpp"

inline void bake(std::string name, bool include_fresnel_dielectrics) {
	rt::SpecularAlbedo albedo;
	albedo.build(256, 256, [&](float alpha, float cosTheta) {
		using namespace rt;
		glm::vec3 wo = glm::vec3(std::sqrt(1.0f - cosTheta * cosTheta), 0.0f, cosTheta);
		glm::vec3 Ng(0.0f, 0.0f, 1.0f);

		//int SampleCount = 100000;
		int SampleCount = 300000;
		Xor64 random;

		OnlineMean<float> mean;
		for (int i = 0; i < SampleCount; ++i) {
			glm::vec3 wi = BeckmannImportanceSampler::sample(&random, alpha, wo, Ng);
			float pdf_omega = BeckmannImportanceSampler::pdf(wi, alpha, wo, Ng);

			glm::vec3 h = glm::normalize(wi + wo);
			float d = D_Beckmann(Ng, h, alpha);
			float g = G2_height_correlated_beckmann(wi, wo, h, Ng, alpha);

			float cos_term_wo = glm::dot(Ng, wo);
			float cos_term_wi = glm::dot(Ng, wi);

			float brdf = d * g / (4.0f * cos_term_wo * cos_term_wi);
			if (include_fresnel_dielectrics) {
				float cosThetaFresnel = glm::dot(h, wo); // == glm::dot(h, wi)
				float f = fresnel_dielectrics(cosThetaFresnel);
				// float f = fresnel_shlick(0.04f, cosThetaFresnel);
				brdf *= f;
			}

			float value = brdf * cos_term_wi / pdf_omega;

			// wiが裏面
			if (glm::dot(wi, Ng) <= 0.0) {
				value = 0.0f;
			}

			if (glm::isfinite(value) == false) {
				value = 0.0f;
			}

			mean.addSample(value);
		}
		return mean.mean();
	});
	albedo.save(ofToDataPath(name + ".xml").c_str());

	// preview
	ofFloatImage image;
	image.allocate(albedo.alphaSize(), albedo.cosThetaSize(), OF_IMAGE_GRAYSCALE);
	float *p = image.getPixels().getPixels();
	for (int j = 0; j < albedo.cosThetaSize(); ++j) {
		for (int i = 0; i < albedo.alphaSize(); ++i) {
			float value = albedo.get(i, j);
			p[j * albedo.alphaSize() + i] = value;
		}
	}
	image.save(name + ".exr");
}

// 両方拡張子を含む
void bake_avg(const char *albedoFile, const char *dstName) {
	rt::SpecularAlbedo albedo;
	albedo.load(ofToDataPath(albedoFile).c_str());
	
	rt::SpecularAvgAlbedo avg;
	avg.build(256, [&](double alpha) {
		return 2.0 * rt::composite_simpson<double>([&](double theta) {
			double cosTheta = cos(theta);
			return albedo.sample(alpha, cosTheta) * cosTheta * sin(theta);
		}, 128, 0.0, glm::pi<double>() * 0.5);
	});
	avg.save(ofToDataPath(dstName).c_str());
}

//--------------------------------------------------------------
void ofApp::setup() {
	// rt::Stopwatch sw;
	// bake("albedo_specular_conductor", false);
	// bake("albedo_specular_dielectrics", true);
	// printf("done %f seconds\n", sw.elapsed());

	// bake_avg("albedo_specular_conductor.xml", "albedo_specular_conductor_avg.xml");
	// bake_avg("albedo_specular_dielectrics.xml", "albedo_specular_dielectrics_avg.xml");
}

//--------------------------------------------------------------
void ofApp::update(){

}

//--------------------------------------------------------------
void ofApp::draw(){

}

//--------------------------------------------------------------
void ofApp::keyPressed(int key){

}

//--------------------------------------------------------------
void ofApp::keyReleased(int key){

}

//--------------------------------------------------------------
void ofApp::mouseMoved(int x, int y ){

}

//--------------------------------------------------------------
void ofApp::mouseDragged(int x, int y, int button){

}

//--------------------------------------------------------------
void ofApp::mousePressed(int x, int y, int button){

}

//--------------------------------------------------------------
void ofApp::mouseReleased(int x, int y, int button){

}

//--------------------------------------------------------------
void ofApp::mouseEntered(int x, int y){

}

//--------------------------------------------------------------
void ofApp::mouseExited(int x, int y){

}

//--------------------------------------------------------------
void ofApp::windowResized(int w, int h){

}

//--------------------------------------------------------------
void ofApp::gotMessage(ofMessage msg){

}

//--------------------------------------------------------------
void ofApp::dragEvent(ofDragInfo dragInfo){ 

}
