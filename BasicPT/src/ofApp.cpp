#include "alembic_loader.hpp"
#include "integrator.hpp"

#include <random>
#include <xmmintrin.h>
#include <pmmintrin.h>
#include <tbb/tbb.h>

#include "ofApp.h"
#include "ofxImGuiLite.hpp"

inline ofPixels toOf(const rt::Image &image) {
	ofPixels pixels;
	pixels.allocate(image.width(), image.height(), OF_IMAGE_COLOR);
	uint8_t *dst = pixels.getPixels();

	double scale = 1.0;
	for (int y = 0; y < image.height(); ++y) {
		for (int x = 0; x < image.width(); ++x) {
			int index = y * image.width() + x;
			const auto &px = *image.pixel(x, y);
			auto L = px.color / (double)px.sample;
			dst[index * 3 + 0] = (uint8_t)glm::clamp(glm::pow(L.x * scale, 1.0 / 2.2) * 255.0, 0.0, 255.99999);
			dst[index * 3 + 1] = (uint8_t)glm::clamp(glm::pow(L.y * scale, 1.0 / 2.2) * 255.0, 0.0, 255.99999);
			dst[index * 3 + 2] = (uint8_t)glm::clamp(glm::pow(L.z * scale, 1.0 / 2.2) * 255.0, 0.0, 255.99999);
		}
	}
	return pixels;
}
inline ofFloatPixels toOfLinear(const rt::Image &image) {
	ofFloatPixels pixels;
	pixels.allocate(image.width(), image.height(), OF_IMAGE_COLOR);
	float *dst = pixels.getPixels();

	for (int y = 0; y < image.height(); ++y) {
		for (int x = 0; x < image.width(); ++x) {
			int index = y * image.width() + x;
			const auto &px = *image.pixel(x, y);
			auto L = px.color / (double)px.sample;
			dst[index * 3 + 0] = L[0];
			dst[index * 3 + 1] = L[1];
			dst[index * 3 + 2] = L[2];
		}
	}
	return pixels;
}

std::shared_ptr<rt::Scene> scene;
std::shared_ptr<rt::PTRenderer> renderer;

bool isPowerOfTwo(uint32_t value)
{
	return value && !(value & (value - 1));
}

//--------------------------------------------------------------
void ofApp::setup() {
	_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
	_MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);

	rt::microcylinder_validation();

	ofxImGuiLite::initialize();

	ofSetVerticalSync(false);

	_camera.setNearClip(0.1);
	_camera.setFarClip(100.0);
	_camera.setDistance(5.0);

	rt::CoupledBRDFConductor::load(
		ofToDataPath("baked/albedo_specular_conductor.bin").c_str(),
		ofToDataPath("baked/albedo_specular_conductor_avg.bin").c_str());
	rt::CoupledBRDFDielectrics::load(
		ofToDataPath("baked/albedo_specular_dielectrics.bin").c_str(),
		ofToDataPath("baked/albedo_specular_dielectrics_avg.bin").c_str());

	rt::CoupledBRDFVelvet::load(
		ofToDataPath("baked/albedo_velvet.bin").c_str(),
		ofToDataPath("baked/albedo_velvet_avg.bin").c_str());

	loadScene();
}
void ofApp::exit() {
	// ofxImGuiLite::shutdown();
}

void ofApp::loadScene() {
	rt::Stopwatch sw;

	scene = std::shared_ptr<rt::Scene>(new rt::Scene());
	rt::loadFromABC(ofToDataPath("cornelbox.abc").c_str(), *scene);
	printf("load scene %f seconds\n", sw.elapsed());

	renderer = std::shared_ptr<rt::PTRenderer>(new rt::PTRenderer(scene));
}
//--------------------------------------------------------------
void ofApp::update() {
}

//--------------------------------------------------------------
void ofApp::draw() {
	ofEnableDepthTest();

	ofClear(0);
	_camera.begin();
	ofPushMatrix();
	ofRotateZ(90.0);
	ofSetColor(255);
	ofDrawGridPlane(1.0);
	ofPopMatrix();

	ofPushMatrix();
	ofDrawAxis(50);
	ofPopMatrix();

	ofSetColor(255);

	if (_showWireframe) {
		for (int i = 0; i < scene->geometries.size(); ++i) {
			static ofMesh mesh;
			mesh.clear();

			auto geometry = scene->geometries[i];
			for (int j = 0; j < geometry.points.size(); ++j) {
				auto p = geometry.points[j].P;
				mesh.addVertex(ofVec3f(p.x, p.y, p.z));
			}
			for (int j = 0; j < geometry.primitives.size(); ++j) {
				auto prim = geometry.primitives[j];
				mesh.addIndex(prim.indices[0]);
				mesh.addIndex(prim.indices[1]);
				mesh.addIndex(prim.indices[2]);
			}
			mesh.drawWireframe();
		}
		{
			auto camera = scene->camera;
			ofSetColor(255);
			auto origin = camera.origin();
			ofDrawSphere(origin.x, origin.y, origin.z, 0.05);

			ofSetColor(0, 0, 255);
			auto frontP = origin + camera.front() * 0.2;
			ofDrawLine(origin.x, origin.y, origin.z, frontP.x, frontP.y, frontP.z);

			ofSetColor(255, 0, 0);
			auto upP = origin + camera.up() * 0.2;
			ofDrawLine(origin.x, origin.y, origin.z, upP.x, upP.y, upP.z);
		}

		{
			rt::Xor64 random;
			int x = ofMap(ofGetMouseX(), 0, ofGetWidth(), 0, scene->camera.imageWidth());
			int y = ofMap(ofGetMouseY(), 0, ofGetHeight(), 0, scene->camera.imageHeight());
			glm::dvec3 o;
			glm::dvec3 d;
			scene->camera.sampleRay(&random, x, y, &o, &d);

			rt::Material m;
			float tmin = 0.0f;
			if (renderer->sceneInterface().intersect(o, d, &m, &tmin)) {
				ofSetColor(255, 0, 0);
				auto p = o + d * (double)tmin;
				ofDrawLine(o.x, o.y, o.z, p.x, p.y, p.z);

				auto pn = p + m->Ng * 0.1;
				ofDrawLine(p.x, p.y, p.z, pn.x, pn.y, pn.z);
			}
			else {
				ofSetColor(255);
				auto p = o + d * 10.0;
				ofDrawLine(o.x, o.y, o.z, p.x, p.y, p.z);
			}
		}
	}

	if(_render) {
		renderer->step();

		ofDisableArbTex();

		if (ofGetFrameNum() % 5 == 0) {
			_image.setFromPixels(toOf(renderer->_image));
		}
		uint32_t n = renderer->stepCount();

		if (32 <= n && isPowerOfTwo(n)) {
			_image.setFromPixels(toOf(renderer->_image));
			char name[64];
			sprintf(name, "%dspp.png", n);
			_image.save(name);
			printf("elapsed %fs\n", ofGetElapsedTimef());
		}

		ofEnableArbTex();
	}

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
	ImGui::SetNextWindowSize(ImVec2(700, 600), ImGuiSetCond_Appearing);
	ImGui::SetNextWindowCollapsed(false, ImGuiSetCond_Appearing);
	ImGui::SetNextWindowBgAlpha(0.5f);

	ImGui::Begin("settings", nullptr);
	ImGui::Checkbox("render", &_render);
	ImGui::Checkbox("show wireframe", &_showWireframe);
	
	ImGui::Text("%d sample, fps = %.3f", renderer->stepCount(), ofGetFrameRate());
	ImGui::Text("%d bad sample nan", renderer->badSampleNanCount());
	ImGui::Text("%d bad sample inf", renderer->badSampleInfCount());
	ImGui::Text("%d bad sample neg", renderer->badSampleNegativeCount());
	ImGui::Text("%d bad sample firefly", renderer->badSampleFireflyCount());
	if (_image.isAllocated()) {
		ofxImGuiLite::image(_image);
	}
	ImGui::End();
}

//--------------------------------------------------------------
void ofApp::keyPressed(int key) {
	if (key == 'f') {
		auto camera = scene->camera;
		ofSetWindowShape(camera.imageWidth(), camera.imageHeight());
		_camera.setNearClip(0.1);
		_camera.setFarClip(100.0);
		_camera.setFov(glm::degrees(camera.setting().fovy));
		_camera.setPosition(camera.origin().x, camera.origin().y, camera.origin().z);

		auto lookAt = camera.origin() + camera.front();
		_camera.lookAt(ofVec3f(lookAt.x, lookAt.y, lookAt.z), ofVec3f(camera.up().x, camera.up().y, camera.up().z));
	}

	if (key == 's') {
		ofFloatImage image = toOfLinear(renderer->_image);
		image.save("pt.exr");
	}

	if (key == 'r') {
		loadScene();
	}
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
