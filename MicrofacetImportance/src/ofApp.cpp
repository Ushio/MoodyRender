#include "ofApp.h"
#include "peseudo_random.hpp"
#include "microfacet.hpp"
//--------------------------------------------------------------
void ofApp::setup(){
	_camera.setNearClip(0.1f);
	_camera.setFarClip(100.0f);
	_camera.setDistance(5.0f);
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
	ofRotateZ(90.0f);
	ofSetColor(128);
	ofDrawGridPlane(1.0f);
	ofPopMatrix();

	ofPushMatrix();
	ofDrawAxis(50);
	ofPopMatrix();

	using namespace rt;

	Xor64 random;
	float alpha = 0.9f;
	float cosTheta = ofMap(ofGetMouseX(), 0, ofGetWidth(), 0, 1);
	glm::vec3 wo = glm::vec3(std::sqrt(1.0f - cosTheta * cosTheta), 0.0f, cosTheta);
	glm::vec3 Ng(0.0f, 0.0f, 1.0f);

	ofSetColor(255, 0, 0);
	ofDrawLine(glm::vec3(), glm::vec3(wo.x, wo.z, wo.y));

	ofMesh mesh;
	mesh.setMode(OF_PRIMITIVE_POINTS);
	for (int i = 0; i < 3000; ++i) {
		glm::vec3 wi = BeckmannImportanceSampler::sample(&random, alpha, wo, Ng);
		mesh.addVertex(glm::vec3(wi.x, wi.z, wi.y));
	}
	ofSetColor(255);
	mesh.draw();

	// 裏側の数値がとても不安定..
	ofPolyline line;
	int N = 3000;
	for (int i = 0; i < N; ++i) {
		float theta = ofMap(i, 0, N - 1, 0, glm::two_pi<float>());
		float cosTheta = cos(theta);
		float sinTheta = sin(theta);
		glm::vec3 wi = glm::vec3(sinTheta, 0.0f, cosTheta);

		double p = BeckmannImportanceSampler::pdf(wi, alpha, wo, Ng);
		if (p < 0) {
			printf("w");
		}
		BeckmannImportanceSampler::pdf(wi, alpha, wo, Ng);

		line.addVertex(glm::vec3(wi.x, wi.z, wi.y) * p);
	}
	line.draw();

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
