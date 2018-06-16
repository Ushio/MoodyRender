#include "ofApp.h"
#include "peseudo_random.hpp"
#include "microfacet.hpp"
#include "material.hpp"

//--------------------------------------------------------------
void ofApp::setup(){
	rt::CoupledBRDFConductor::load(ofToDataPath("baked/albedo_specular_conductor.xml").c_str(), ofToDataPath("baked/albedo_specular_conductor_avg.xml").c_str());
	rt::CoupledBRDFDielectrics::load(ofToDataPath("baked/albedo_specular_dielectrics.xml").c_str(), ofToDataPath("baked/albedo_specular_dielectrics_avg.xml").c_str());

	_camera.setNearClip(0.1);
	_camera.setFarClip(100.0);
	_camera.setDistance(5.0);
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
	ofSetColor(128);
	ofDrawGridPlane(1.0);
	ofPopMatrix();

	ofPushMatrix();
	ofDrawAxis(50);
	ofPopMatrix();

	using namespace rt;

	Xor64 random;
	double alpha = 0.3;
	double cosTheta = ofMap(ofGetMouseX(), 0, ofGetWidth(), 0, 1);
	glm::dvec3 wo = glm::dvec3(std::sqrt(1.0 - cosTheta * cosTheta), 0.0, cosTheta);
	glm::dvec3 Ng(0.0, 0.0, 1.0);

	ofSetColor(255, 0, 0);
	ofDrawLine(glm::vec3(), glm::vec3(wo.x, wo.z, wo.y));

	ofMesh mesh;
	mesh.setMode(OF_PRIMITIVE_POINTS);
	for (int i = 0; i < 3000; ++i) {
		double theta = CoupledBRDFConductor::sampler().sampleTheta(alpha, &random);
		// double theta = random.uniform(0.0, glm::pi<double>() * 0.5);
		glm::dvec3 sample = polar_to_cartesian(theta, random.uniform(0.0, glm::two_pi<double>()));
		ArbitraryBRDFSpace space(Ng);
		glm::dvec3 wi = space.localToGlobal(sample);

		// glm::dvec3 wi = BeckmannImportanceSampler::sample(&random, alpha, wo, Ng);
		mesh.addVertex(glm::dvec3(wi.x, wi.z, wi.y));
	}
	ofSetColor(255);
	mesh.draw();

	//// 裏側の数値がとても不安定..
	//ofPolyline line;
	//int N = 3000;
	//for (int i = 0; i < N; ++i) {
	//	double theta = ofMap(i, 0, N - 1, 0, glm::two_pi<double>());
	//	double cosTheta = cos(theta);
	//	double sinTheta = sin(theta);
	//	glm::dvec3 wi = glm::dvec3(sinTheta, 0.0, cosTheta);

	//	double p = BeckmannImportanceSampler::pdf(wi, alpha, wo, Ng);
	//	if (p < 0) {
	//		printf("w");
	//	}
	//	BeckmannImportanceSampler::pdf(wi, alpha, wo, Ng);

	//	line.addVertex(glm::dvec3(wi.x, wi.z, wi.y) * p);
	//}
	//line.draw();

	//_camera.end();

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
