#include "ofApp.h"
#include "bicubic.hpp"

#include <algorithm>

// 基本的な擬似乱数
namespace rt {
	struct PeseudoRandom {
		virtual ~PeseudoRandom() {}

		virtual uint32_t generate() = 0;
		virtual double uniform() = 0;
		virtual double uniform(double a, double b) = 0;

		virtual float uniformf() { return (float)uniform(); }
		virtual float uniformf(float a, float b) { return (float)(uniform(a, b)); }
	};
	struct Xor : public PeseudoRandom {
		Xor() {

		}
		Xor(uint32_t seed) {
			_y = std::max(seed, 1u);
		}

		// 0 <= x <= 0x7FFFFFFF
		uint32_t generate() {
			_y = _y ^ (_y << 13); _y = _y ^ (_y >> 17);
			uint32_t value = _y = _y ^ (_y << 5); // 1 ~ 0xFFFFFFFF(4294967295
			return value >> 1;
		}
		// 0.0 <= x < 1.0
		double uniform() {
			return double(generate()) / double(0x80000000);
		}
		double uniform(double a, double b) {
			return a + (b - a) * double(uniform());
		}
	public:
		uint32_t _y = 2463534242;
	};
}

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
void ofApp::draw() {
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

#if 0
	rt::Xor random;
	float image[4];
	for (int i = 0; i < 4; ++i) {
		image[i] = random.uniformf();
	}

	ofSetColor(255);
	for (int i = 0; i < 4; ++i) {
		ofDrawSphere(i, image[i], 0.05f);
	}

	ofMesh mesh;
	mesh.setMode(OF_PRIMITIVE_POINTS);

	int N = 1000; 
	for (int i = 0; i < N; ++i) {
		float x = ofMap(i, 0, N - 1, -0.2f, 1.2f);
		float value = bicubic_1d(x, 4, [&](int x) { return image[x]; });
		mesh.addVertex(glm::vec3(x * 3, value, 0.0f));
	}
	mesh.draw();
#else
	const int X_SIZE = 6;
	const int Y_SIZE = 6;

	rt::Xor random;
	float image[Y_SIZE][X_SIZE];
	for (int z = 0; z < Y_SIZE; ++z) {
		for (int x = 0; x < X_SIZE; ++x) {
			image[z][x] = random.uniformf();
		}
	}

	ofSetColor(255);
	for (int z = 0; z < Y_SIZE; ++z) {
		for (int x = 0; x < X_SIZE; ++x) {
			ofDrawSphere(x, image[z][x], z, 0.05f);
		}
	}

	ofMesh mesh;
	mesh.setMode(OF_PRIMITIVE_POINTS);

	int N = 500;
	for (int yi = 0; yi < N; ++yi) {
		float y = ofMap(yi, 0, N - 1, -0.5f, 1.5f);

		for (int xi = 0; xi < N; ++xi) {
			float x = ofMap(xi, 0, N - 1, -0.5f, 1.5f);

			float value = bicubic_2d(x, y, X_SIZE, Y_SIZE, [&](int x, int y) { return image[y][x]; });
			mesh.addVertex(glm::vec3(x * (X_SIZE - 1), value, y * (Y_SIZE - 1)));
		}
	}
	mesh.draw();
#endif

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
