#include "ofApp.h"

#include "serializable_buffer.hpp"
#include "peseudo_random.hpp"
#include "microfacet.hpp"
#include "online.hpp"
#include "stopwatch.hpp"
#include "composite_simpson.hpp"
#include "simpson_helper.hpp"
#include "linear_transform.hpp"

static const int kBakeResolution = 64;

inline void bake(std::string name, bool include_fresnel_dielectrics) {
	rt::SpecularAlbedo albedo;
	albedo.build(kBakeResolution, kBakeResolution, [&](double alpha, double cosTheta) {
		using namespace rt;
		glm::dvec3 wo = glm::dvec3(std::sqrt(1.0 - cosTheta * cosTheta), 0.0, cosTheta);
		glm::dvec3 Ng(0.0, 0.0, 1.0);

		//int SampleCount = 100000;
		int SampleCount = 300000;
		Xor64 random;

		OnlineMean<double> mean;
		for (int i = 0; i < SampleCount; ++i) {
			glm::dvec3 wi = VCavityBeckmannVisibleNormalSampler::sample(&random, alpha, wo, Ng);
			double pdf_omega = VCavityBeckmannVisibleNormalSampler::pdf(wi, alpha, wo, Ng);

			glm::dvec3 h = glm::normalize(wi + wo);
			double d = D_Beckmann(Ng, h, alpha);
			double g = G2_v_cavity(wo, wi, h, Ng);

			double cos_term_wo = glm::dot(Ng, wo);
			double cos_term_wi = glm::dot(Ng, wi);

			double brdf = d * g / (4.0 * cos_term_wo * cos_term_wi);
			if (include_fresnel_dielectrics) {
				double cosThetaFresnel = glm::dot(h, wo); // == glm::dot(h, wi)
				double f = fresnel_dielectrics(cosThetaFresnel);
				brdf *= f;
			}

			double value = brdf * cos_term_wi / pdf_omega;

			// wiが裏面
			if (glm::dot(wi, Ng) <= 0.0) {
				value = 0.0;
			}

			if (glm::isfinite(value) == false) {
				value = 0.0;
			}

			mean.addSample(value);
		}
		return mean.mean();
	});
	saveAsBinary(albedo, ofToDataPath(name + ".bin").c_str());

	// preview
	ofFloatImage image;
	image.allocate(albedo.alphaSize(), albedo.cosThetaSize(), OF_IMAGE_GRAYSCALE);
	float *p = image.getPixels().getPixels();
	for (int j = 0; j < albedo.cosThetaSize(); ++j) {
		for (int i = 0; i < albedo.alphaSize(); ++i) {
			double value = albedo.get(i, j);
			p[j * albedo.alphaSize() + i] = value;
		}
	}
	image.save(name + ".exr");
}

// 両方拡張子を含む
void bake_avg(const char *albedoFile, const char *dstName) {
	rt::SpecularAlbedo albedo;
	loadFromBinary(albedo, ofToDataPath(albedoFile).c_str());
	
	rt::SpecularAvgAlbedo avg;
	avg.build(kBakeResolution, [&](double alpha) {
		return 2.0 * rt::composite_simpson<double>([&](double theta) {
			double cosTheta = cos(theta);
			return albedo.sample(alpha, cosTheta) * cosTheta * sin(theta);
		}, kBakeResolution, 0.0, glm::pi<double>() * 0.5);
	});
	saveAsBinary(avg, ofToDataPath(dstName).c_str());
}


//--------------------------------------------------------------
void ofApp::setup() {
	 //rt::Stopwatch sw;
	 //bake("albedo_specular_conductor", false);
	 //bake("albedo_specular_dielectrics", true);
	 //printf("done %f seconds\n", sw.elapsed());

	 bake_avg("albedo_specular_conductor.bin", "albedo_specular_conductor_avg.bin");
	 // bake_avg("albedo_specular_dielectrics.bin", "albedo_specular_dielectrics_avg.bin");

	ofSetVerticalSync(false);

	_camera.setNearClip(0.1);
	_camera.setFarClip(100.0);
	_camera.setDistance(5.0);

	// rt::CoupledBRDFConductor::load(ofToDataPath("albedo_specular_conductor.xml").c_str(), ofToDataPath("albedo_specular_conductor_avg.xml").c_str());
	// rt::CoupledBRDFDielectrics::load(ofToDataPath("albedo_specular_dielectrics.xml").c_str(), ofToDataPath("albedo_specular_dielectrics_avg.xml").c_str());
}

//--------------------------------------------------------------
void ofApp::update(){

}

//--------------------------------------------------------------
void ofApp::draw(){
	ofEnableDepthTest();

	ofClear(0);
	_camera.begin();
	ofPushMatrix();
	ofRotateZ(90.0);
	ofSetColor(64);
	ofDrawGridPlane(1.0);
	ofPopMatrix();

	ofPushMatrix();
	ofDrawAxis(50);
	ofPopMatrix();

	

	//{
	//	ofMesh mesh;
	//	mesh.setMode(OF_PRIMITIVE_POINTS);

	//	int N = 256;
	//	std::vector<double> values(N);

	//	rt::LinearTransform<double> transform(0, N - 1, 0, 1);
	//	for (int xi = 0; xi < N; ++xi) {
	//		double alpha = transform(xi);

	//		for (int zi = 0; zi < N; ++zi) {
	//			double cosTheta = transform(zi);
	//			// CoupledBRDFConductor
	//			// CoupledBRDFDielectrics

	//			//double value = (1.0 - rt::CoupledBRDFDielectrics::specularAlbedo().sample(alpha, cosTheta)) * cosTheta;
	//			//mesh.addVertex(glm::dvec3(alpha, value, cosTheta));

	//			values[zi] = (1.0 - rt::CoupledBRDFConductor::specularAlbedo().sample(alpha, cosTheta)) * cosTheta;
	//			
	//			// values[zi] = (1.0 - rt::CoupledBRDFDielectrics::specularAlbedo().sample(alpha, cosTheta)) * cosTheta;
	//		}

	//		double maxValue = *std::max_element(values.begin(), values.end());

	//		for (int zi = 0; zi < N; ++zi) {
	//			double cosTheta = transform(zi);
	//			mesh.addVertex(glm::dvec3(alpha, values[zi], cosTheta));
	//		}
	//	}

	//	ofSetColor(255);
	//	mesh.draw();
	//}

	_camera.end();

	ofDisableDepthTest();
	ofSetColor(255);

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
