#include "ofApp.h"

#include <glm/glm.hpp>
#include <glm/ext.hpp>
#include <functional>

#include "peseudo_random.hpp"

inline float bicubic_kernel(float t, float tt, float ttt, float f0, float f1, float f2, float f3) {
	glm::mat4 m = {
		0.0f, -1.0f, 2.0f, -1.0f,
		2.0f, 0.0f, -5.0f, 3.0f,
		0.0f, 1.0f, 4.0f, -3.0f,
		0.0f, 0.0f, -1.0f, 1.0f
	};
	return 0.5f * glm::dot(glm::vec4(1.0f, t, tt, ttt), m * glm::vec4(f0, f1, f2, f3));
}
inline float bicubic_kernel(float t, float f0, float f1, float f2, float f3) {
	float tt = t * t;
	float ttt = tt * t;
	return bicubic_kernel(t, tt, ttt, f0, f1, f2, f3);
}

// f(0.0) = sample(0)
// f(1.0) = sample(size - 1)
inline float bicubic_1d(float x /* 0.0 => 1.0 */, int size, std::function<float(int)> sample) {
	float index_f = x * (size - 1);
	int index1 = (int)std::floor(index_f);

	float values[4];
	for (int i = 0; i < 4; ++i) {
		int index = index1 - 1 + i;
		index = glm::clamp(index, 0, size - 1);
		values[i] = sample(index);
	}

	float t = index_f - index1;
	return bicubic_kernel(t, values[0], values[1], values[2], values[3]);
}

inline float bicubic_2d(float x /* 0.0 => 1.0 */, float y /* 0.0 => 1.0 */, int sizex, int sizey, std::function<float(int, int)> sample) {
	float index_xf = x * (sizex - 1);
	float index_yf = y * (sizey - 1);
	int index1x = (int)std::floor(index_xf);
	int index1y = (int)std::floor(index_yf);
	float tx = index_xf - index1x;
	float ty = index_yf - index1y;

	float values[4][4];
	for (int y = 0; y < 4; ++y) {
		int yi = index1y - 1 + y;
		yi = glm::clamp(yi, 0, sizey - 1);
		for (int x = 0; x < 4; ++x) {
			int xi = index1x - 1 + x;
			xi = glm::clamp(xi, 0, sizex - 1);
			values[y][x] = sample(xi, yi);
		}
	}

	float ttx = tx * tx;
	float tttx = ttx * tx;
	return bicubic_kernel(ty, 
		bicubic_kernel(tx, ttx, tttx, values[0][0], values[0][1], values[0][2], values[0][3]),
		bicubic_kernel(tx, ttx, tttx, values[1][0], values[1][1], values[1][2], values[1][3]),
		bicubic_kernel(tx, ttx, tttx, values[2][0], values[2][1], values[2][2], values[2][3]),
		bicubic_kernel(tx, ttx, tttx, values[3][0], values[3][1], values[3][2], values[3][3])
	);
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
	const int XDIM = 6;
	const int YDIM = 6;

	rt::Xor random;
	float image[YDIM][XDIM];
	for (int z = 0; z < YDIM; ++z) {
		for (int x = 0; x < XDIM; ++x) {
			image[z][x] = random.uniformf();
		}
	}

	ofSetColor(255);
	for (int z = 0; z < YDIM; ++z) {
		for (int x = 0; x < XDIM; ++x) {
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

			float value = bicubic_2d(x, y, XDIM, YDIM, [&](int x, int y) { return image[y][x]; });
			mesh.addVertex(glm::vec3(x * (XDIM - 1), value, y * (YDIM - 1)));
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
