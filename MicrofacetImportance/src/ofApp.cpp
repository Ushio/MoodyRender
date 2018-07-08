#include "ofApp.h"
#include "ofxImGuiLite.hpp"

#include "peseudo_random.hpp"
#include "microfacet.hpp"
#include "material.hpp"
#include "simpson_helper.hpp"


inline double G_dist_1_smith(const glm::dvec3 &omega, const glm::dvec3 &Ng, double alpha) {
	double cosThetaO = glm::abs(glm::dot(Ng, omega));
	double tanThetaO = std::sqrt(1.0 - cosThetaO * cosThetaO) / cosThetaO;
	double a = 1.0 / (alpha * tanThetaO);
	auto lambda = (std::erf(a) - 1.0) * 0.5 + std::exp(-a * a) / (2.0 * a * std::sqrt(glm::pi<double>()));
	return 1.0 / (1.0 + lambda);
}
inline double G_1_smith(const glm::dvec3 &omega, const glm::dvec3 &omega_m, const glm::dvec3 &Ng, double alpha) {
	return rt::chi_plus(glm::dot(omega, omega_m)) * G_dist_1_smith(omega, Ng, alpha);
}

inline double G1_v_cavity_basic(glm::dvec3 omega, glm::dvec3 H, glm::dvec3 N) {
	double a = 2.0 * glm::dot(N, H) * glm::dot(N, omega) / glm::max(glm::dot(omega, H), 0.0);
	return glm::min(a, 1.0);
}
inline double G1_v_cavity_ext(glm::dvec3 omega, glm::dvec3 H, glm::dvec3 N) {
	double chi = rt::chi_plus(glm::dot(omega, H) / glm::dot(N, omega));
	double a = 2.0 * glm::abs(glm::dot(N, H)) * glm::abs(glm::dot(N, omega)) / glm::abs(glm::dot(omega, H));
	return chi * glm::min(a, 1.0);
}
//--------------------------------------------------------------
void ofApp::setup(){
	ofxImGuiLite::initialize();

	using namespace rt;

	// visible normalが、
	// glm::max(glm::dot(wo, wm), 0.0) でクランプされていたからうまくいっていただけだ。

	//rt::Xor64 random;
	//glm::dvec3 Ng(0, 0, 1);

	//for (int j = 0; j < 32; ++j) {
	//	// alphaが小さい場合、simpsonによる積分が適さない
	//	double alpha = random.uniform(0.1, 1.0);
	//	glm::dvec3 wo = LambertianSampler::sample(&random, Ng);
	//	double cosThetaO = glm::dot(wo, Ng);

	//	double integral = hemisphere_composite_simpson<double>([&](double theta, double phi) {
	//		glm::dvec3 wm = rt::polar_to_cartesian((double)theta, (double)phi);
	//		// double g1 = G1_v_cavity(wo, wm, Ng);
	//		// double g1 = G_1_smith(wo, wm, Ng, alpha);
	//		double g1 = G1_v_cavity_ext(wo, wm, Ng);
	//		double value = g1 * glm::max(glm::dot(wo, wm), 0.0) * D_Beckmann(Ng, wm, alpha) / cosThetaO;
	//		return value;
	//	}, 500);
	//	
	//	printf("[%d] alpha = %.5f / cosThetaO = %.5f  / %.5f\n", j, alpha, cosThetaO, integral);
	//}

	rt::CoupledBRDFConductor::load(ofToDataPath("baked/albedo_specular_conductor.bin").c_str(), ofToDataPath("baked/albedo_specular_conductor_avg.bin").c_str());
	rt::CoupledBRDFDielectrics::load(ofToDataPath("baked/albedo_specular_dielectrics.bin").c_str(), ofToDataPath("baked/albedo_specular_dielectrics_avg.bin").c_str());

	_camera.setNearClip(0.1);
	_camera.setFarClip(100.0);
	_camera.setDistance(5.0);
}

//--------------------------------------------------------------
void ofApp::update(){

}

//--------------------------------------------------------------
void ofApp::draw(){
	using namespace rt;

	static float theta = 0.0;
	static float alpha = 0.4;
	static bool isVisibleNormal = true;

	ofEnableDepthTest();

	ofClear(0);
	_camera.begin();
	ofPushMatrix();
	ofRotateYDeg(90.0);
	// ofRotateYDeg(90.0);
	ofSetColor(128);
	ofDrawGridPlane(1.0);
	ofPopMatrix();

	ofPushMatrix();
	ofDrawAxis(50);
	ofPopMatrix();


	Xor64 random;

	double cosTheta = cos(theta);
	double sinTheta = sin(theta);
	glm::dvec3 wo = glm::dvec3(sinTheta, 0.0, cosTheta);
	glm::dvec3 Ng(0.0, 0.0, 1.0);

	ofSetColor(255, 0, 0);
	ofDrawLine(glm::vec3(), glm::vec3(wo.x, wo.z, wo.y));

	{
		ofMesh mesh;
		mesh.setMode(OF_PRIMITIVE_POINTS);
		for (int i = 0; i < 3000; ++i) {
			glm::dvec3 wi = VelvetSampler::sample(&random, alpha, wo, Ng);
			mesh.addVertex(glm::dvec3(wi.x, wi.z, wi.y));
		}
		ofSetColor(255);
		mesh.draw();

	}
	{
		ofPolyline line;
		int N = 3000;
		for (int i = 0; i < N; ++i) {
			double theta = ofMap(i, 0, N - 1, 0, glm::two_pi<double>());
			double cosTheta = cos(theta);
			double sinTheta = sin(theta);
			glm::dvec3 wi = glm::dvec3(sinTheta, 0.0, cosTheta);
			// glm::dvec3 half = glm::normalize(wi + wo);
			double p = D_velvet(Ng, wi, alpha);

			line.addVertex(glm::dvec3(wi.x, wi.z, wi.y) * p);
		}
		line.draw();
	}

	//ofMesh mesh;
	//mesh.setMode(OF_PRIMITIVE_POINTS);
	//for (int i = 0; i < 3000; ++i) {
	//	double theta = CoupledBRDFConductor::sampler().sampleTheta(alpha, &random);
	//	glm::dvec3 sample = polar_to_cartesian(theta, random.uniform(0.0, glm::two_pi<double>()));
	//	ArbitraryBRDFSpace space(Ng);
	//	glm::dvec3 wi = space.localToGlobal(sample);
	//	mesh.addVertex(glm::dvec3(wi.x, wi.z, wi.y));
	//}
	//ofSetColor(255);
	//mesh.draw();

	//{
	//	ofPolyline line;
	//	int N = 3000;
	//	for (int i = 0; i < N; ++i) {
	//		double theta = ofMap(i, 0, N - 1, 0, glm::two_pi<double>());
	//		double cosTheta = cos(theta);
	//		double sinTheta = sin(theta);
	//		glm::dvec3 wi = glm::dvec3(sinTheta, 0.0, cosTheta);

	//		double p = BeckmannImportanceSampler::pdf(wi, alpha, wo, Ng);
	//		BeckmannImportanceSampler::pdf(wi, alpha, wo, Ng);

	//		line.addVertex(glm::dvec3(wi.x, wi.z, wi.y) * p);
	//	}
	//	line.draw();
	//}


	//ofMesh mesh;
	//mesh.setMode(OF_PRIMITIVE_POINTS);
	//for (int i = 0; i < 3000; ++i) {
	//	//double theta = CoupledBRDFConductor::sampler().sampleTheta(alpha, &random);
	//	//glm::dvec3 sample = polar_to_cartesian(theta, random.uniform(0.0, glm::two_pi<double>()));
	//	//ArbitraryBRDFSpace space(Ng);
	//	//glm::dvec3 wi = space.localToGlobal(sample);

	//	if (isVisibleNormal) {
	//		glm::dvec3 wi = VCavityBeckmannVisibleNormalSampler::sample(&random, alpha, wo, Ng);
	//		mesh.addVertex(glm::dvec3(wi.x, wi.z, wi.y));
	//	}
	//	else {
	//		glm::dvec3 wi = BeckmannImportanceSampler::sample(&random, alpha, wo, Ng);
	//		mesh.addVertex(glm::dvec3(wi.x, wi.z, wi.y));
	//	}
	//}
	//ofSetColor(255);
	//mesh.draw();

	//double cosTheta = cos(theta);
	//double sinTheta = sin(theta);
	//glm::dvec3 wo = glm::dvec3(sinTheta, 0.0, cosTheta);
	//glm::dvec3 Ng(0.0, 0.0, 1.0);
	//ofPolyline line;
	//int N = 3000;
	//for (int i = 0; i < N; ++i) {
	//	double theta = ofMap(i, 0, N - 1, 0, glm::two_pi<double>());
	//	double cosTheta = cos(theta);
	//	double sinTheta = sin(theta);
	//	glm::dvec3 wi = glm::dvec3(sinTheta, 0.0, cosTheta);

	//	double p = BeckmannImportanceSampler::pdf(wi, alpha, wo, Ng);
	//	BeckmannImportanceSampler::pdf(wi, alpha, wo, Ng);

	//	line.addVertex(glm::dvec3(wi.x, wi.z, wi.y) * p);
	//}
	//line.draw();

	//double cosTheta = cos(theta);
	//double sinTheta = sin(theta);
	//glm::dvec3 wo = glm::dvec3(sinTheta, 0.0, cosTheta);
	//glm::dvec3 Ng(0.0, 0.0, 1.0);

	//{
	//	ofMesh mesh;
	//	mesh.setMode(OF_PRIMITIVE_POINTS);
	//	hemisphere_composite_simpson<double>([&](double theta, double phi) {
	//		glm::dvec3 wm = rt::polar_to_cartesian((double)theta, (double)phi);
	//		mesh.addVertex(glm::dvec3(wm.x, wm.z, wm.y));
	//		return 0;
	//	}, 100);

	//	mesh.draw();
	//}

	
	//ofSetColor(255, 0, 0);
	//ofDrawLine(glm::vec3(), glm::dvec3(wo.x, wo.z, wo.y));

	//ofPolyline line;
	//int N = 3000;
	//for (int i = 0; i < N; ++i) {
	//	double theta = ofMap(i, 0, N - 1, 0, glm::two_pi<double>());
	//	// double theta = ofMap(i, 0, N - 1, 0, glm::pi<double>()) - glm::pi<double>() * 0.5;

	//	double cosTheta = cos(theta);
	//	double sinTheta = sin(theta);
	//	glm::dvec3 wm = glm::dvec3(sinTheta, 0.0, cosTheta);

	//	// double g1 = G1_v_cavity(wo, wm, Ng);
	//	// double g1 = G_1_smith(wo, wm, Ng, 0.5);
	//	// double g1 = G1_v_cavity_ext(wo, wm, Ng);
	//	double g1 = chi_plus(glm::dot(wo, wm) / glm::dot(wo, Ng));

	//	line.addVertex(glm::dvec3(wm.x, wm.z, wm.y) * g1);
	//}

	//ofSetColor(255);
	//line.draw();

	_camera.end();

	ofDisableDepthTest();
	ofSetColor(255);

	ofxImGuiLite::ScopedImGui imgui;

	// camera control                                          for control clicked problem
	if (ImGui::IsWindowHovered(ImGuiHoveredFlags_AnyWindow) || (ImGui::IsAnyWindowFocused() && ImGui::IsAnyMouseDown())) {
		_camera.disableMouseInput();
	}
	else {
		_camera.enableMouseInput();
	}

	ImGui::SetNextWindowPos(ImVec2(20, 20), ImGuiSetCond_Appearing);
	ImGui::SetNextWindowSize(ImVec2(300, 300), ImGuiSetCond_Appearing);
	ImGui::SetNextWindowCollapsed(false, ImGuiSetCond_Appearing);
	ImGui::SetNextWindowBgAlpha(0.5f);

	ImGui::Begin("settings", nullptr);
	ImGui::SliderFloat("alpha", &alpha, 0.0f, 1.0f);
	ImGui::SliderFloat("theta", &theta, 0.0f, glm::pi<float>());

	ImGui::Checkbox("visible normal", &isVisibleNormal);

	ImGui::End();
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
