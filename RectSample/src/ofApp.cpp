#include "ofApp.h"
#include "ofxImGuiLite.hpp"

#include "peseudo_random.hpp"
#include "microfacet.hpp"
#include "material.hpp"
#include "simpson_helper.hpp"

#include "randomsampler.hpp"

#include "online.hpp"

namespace rt {
	struct QuadPlane {
		glm::dvec3 s;
		glm::dvec3 ex;
		glm::dvec3 ey;

		glm::dvec3 sample(double u, double v) const {
			return s + ex * u + ey * v;
		}
		double area() const {
			return glm::length(ex) * glm::length(ey);
		}
		glm::dvec3 normal() const {
			return glm::normalize(glm::cross(ex, ey));
		}
	};

	class SphereAreaSampler {
	public:
		SphereAreaSampler(const QuadPlane &quadPlane, const glm::dvec3 &o) {
			_o = o;

			double exLength = glm::length(quadPlane.ex);
			double eyLength = glm::length(quadPlane.ey);
			_x = quadPlane.ex / exLength;
			_y = quadPlane.ey / eyLength;
			_z = glm::cross(_x, _y);

			glm::dvec3 d = quadPlane.s - o;
			_x0 = glm::dot(d, _x);
			_y0 = glm::dot(d, _y);
			_x1 = _x0 + exLength;
			_y1 = _y0 + eyLength;

			_z0 = glm::dot(d, _z);

			// z flip
			if (_z0 > 0.0) {
				_z0 *= -1.0;
				_z *= -1.0;
			}

			_v[0][0] = glm::vec3(_x0, _y0, _z0);
			_v[0][1] = glm::vec3(_x0, _y1, _z0);
			_v[1][0] = glm::vec3(_x1, _y0, _z0);
			_v[1][1] = glm::vec3(_x1, _y1, _z0);

			auto calculate_n = [](glm::vec3 a, glm::vec3 b) {
				glm::vec3 c = glm::cross(a, b);
				return glm::normalize(c);
			};

			_n[0] = calculate_n(_v[0][0], _v[1][0]);
			_n[1] = calculate_n(_v[1][0], _v[1][1]);
			_n[2] = calculate_n(_v[1][1], _v[0][1]);
			_n[3] = calculate_n(_v[0][1], _v[0][0]);

			_g[0] = std::acos(glm::dot(-_n[0], _n[1]));
			_g[1] = std::acos(glm::dot(-_n[1], _n[2]));
			_g[2] = std::acos(glm::dot(-_n[2], _n[3]));
			_g[3] = std::acos(glm::dot(-_n[3], _n[0]));
		}

		double solidAngle() const {
			return _g[0] + _g[1] + _g[2] + _g[3] - glm::two_pi<double>();
		}

		double minPhiU() const {
			return -_g[2] - _g[3] + glm::two_pi<double>();
		}
		double maxPhiU() const {
			return solidAngle() -_g[2] - _g[3] + glm::two_pi<double>();
		}
		glm::dvec3 sample(double u, double v) const {
			double AQ = solidAngle();

			double phi_u = u * AQ - _g[2] - _g[3] + glm::two_pi<double>();

			auto safeSqrt = [](double x) {
				return std::sqrt(std::max(x, 0.0));
			};

			double b0 = _n[0].z;
			double b1 = _n[2].z;
			double fu = (std::cos(phi_u) * b0 - b1) / std::sin(phi_u);
			double cu = std::copysign(1.0, fu) / safeSqrt(fu * fu + b0 * b0);
			double xu = -cu * _z0 / safeSqrt(1.0 - cu * cu);

			double d = std::sqrt(xu * xu + _z0 * _z0);
			double h0 = _y0 / std::sqrt(d * d + _y0 * _y0);
			double h1 = _y1 / std::sqrt(d * d + _y1 * _y1);
			double hv = glm::mix(h0, h1, v);
			double yv = hv * d / safeSqrt(1.0 - hv * hv);
			return _o + xu * _x + yv * _y + _z0 * _z;
		}

		glm::dvec3 _o;
		glm::dvec3 _x;
		glm::dvec3 _y;
		glm::dvec3 _z;

		// local coord : 
		double _x0;
		double _x1;
		double _y0;
		double _y1;
		double _z0;

		glm::dvec3 _v[2][2];
		glm::dvec3 _n[4];
		double _g[4];
	};
}

inline void drawPlane(rt::QuadPlane plane) {
	// これで同じ現象ということは、xが間違っているらしい
	// plane.ey *= 0.5;

	ofSetColor(255, 0, 0);
	ofDrawLine((glm::vec3)plane.s, plane.s + plane.ex);
	ofSetColor(0, 255, 0);
	ofDrawLine((glm::vec3)plane.s, plane.s + plane.ey);

	rt::SphereAreaSampler sampler(plane, glm::dvec3(0.0));

	ofSetColor(255, 0, 0);
	ofDrawLine((glm::vec3)sampler._o, sampler._o + sampler._x);
	ofSetColor(0, 255, 0);
	ofDrawLine((glm::vec3)sampler._o, sampler._o + sampler._y);
	ofSetColor(0, 0, 255);
	ofDrawLine((glm::vec3)sampler._o, sampler._o + sampler._z);

	static rt::XoroshiroPlus128 random;
	ofMesh mesh;
	mesh.setMode(OF_PRIMITIVE_POINTS);
	for (int i = 0; i < 3000; ++i) {
		double u = random.uniform();
		double v = random.uniform();

		mesh.addVertex(plane.sample(u, v));
		mesh.addVertex(sampler.sample(u, v));
	}
	ofSetColor(255);
	mesh.draw();
}

rt::QuadPlane g_bugPlane;

inline void sampleTest(int seed) {
	rt::XoroshiroPlus128 random(seed);
	glm::dvec3 o = glm::dvec3(random.uniform(-1.0, 1.0), random.uniform(-1.0, 1.0), random.uniform(-1.0, 1.0));

	rt::QuadPlane quadPlane;
	quadPlane.s.x = random.uniform(-1.0, 1.0);
	quadPlane.s.y = random.uniform(3.0, 4.0);
	quadPlane.s.z = random.uniform(-1.0, 1.0);
	
	quadPlane.ex = rt::sample_on_unit_sphere(&random);
	quadPlane.ey = rt::ArbitraryBRDFSpace(quadPlane.ex).xaxis;
	quadPlane.ex *= random.uniform(1.0, 2.0);
	quadPlane.ey *= random.uniform(1.0, 2.0);

	if (seed == 1) {
		g_bugPlane = quadPlane;
	}

	rt::SphereAreaSampler sampler(quadPlane, o);
	double sr = sampler.solidAngle();

	rt::OnlineMean<double> Lo1;
	rt::OnlineMean<double> Lo2;
	glm::dvec3 Ng(0, 1, 0);
	glm::dvec3 LN = quadPlane.normal();
	

	printf("-- (%d) --\n", seed);
	printf("sr: %.4f\n", sr);
	printf("area: %.4f\n", quadPlane.area());

	for (int i = 0; i < 100000; ++i) {
		double Li = 10.0;
		double brdf = 1.0 / glm::pi<double>();
		{
			glm::dvec3 sample = quadPlane.sample(random.uniform(), random.uniform());
			glm::dvec3 wi = glm::normalize(sample - o);

			double G = glm::dot(Ng, wi) * glm::abs(glm::dot(LN, -wi)) / glm::distance2(sample, o);
			double pdf_area = 1.0 / quadPlane.area();
			double value = Li * brdf * G / pdf_area;
			Lo1.addSample(value);
		}

		{
			glm::dvec3 sample = sampler.sample(random.uniform(), random.uniform());
			glm::dvec3 wi = glm::normalize(sample - o);
			double pdf_sr = 1.0 / sampler.solidAngle();
			double value = Li * brdf * glm::dot(Ng, wi) / pdf_sr;
			Lo2.addSample(value);
		}
		if (i % 10000 == 0) {
			// printf("%.5f, %.5f, d = %.5f\n", Lo1.mean(), Lo2.mean(), Lo1.mean() - Lo2.mean());
		}
	}
	printf("%.5f, %.5f, d = %.5f\n", Lo1.mean(), Lo2.mean(), Lo1.mean() - Lo2.mean());
}

//--------------------------------------------------------------
void ofApp::setup(){
	ofxImGuiLite::initialize();

	_camera.setNearClip(0.1);
	_camera.setFarClip(100.0);
	_camera.setDistance(5.0);

	//double inf = INFINITY;
	//printf("%f", sqrt(inf * inf));  // "nan"
	//printf("%f", 3.0f / sqrt(inf));  // "nan"

	// テスト
	for (int i = 0; i < 10; ++i) {
		sampleTest(i);
	}
}

//--------------------------------------------------------------
void ofApp::update(){

}


//--------------------------------------------------------------
void ofApp::draw(){
	using namespace rt;

	ofEnableDepthTest();

	ofClear(0);
	_camera.begin();
	ofPushMatrix();
	// ofRotateYDeg(90.0);
	ofRotateZDeg(90.0);
	ofSetColor(128);
	ofDrawGridPlane(1.0);
	ofPopMatrix();

	ofPushMatrix();
	ofDrawAxis(50);
	ofPopMatrix();

	ofSetCircleResolution(100);



	static float exLen = 2.0f;
	static float eyLen = 1.0f;

	// longitude
	static float alpha = 0.0f;

	// polar
	static float beta = 0.0f;

	static bool bang = false;

	auto fromPolar = [](float a, float b) {
		double sinB = std::sin(b);
		glm::dvec3 v = {
			std::sin(a) * sinB,
			std::cos(b),
			std::cos(a) * sinB
		};
		return v;
	};

	ofSetColor(ofColor::orange);
	ofDrawLine(vec3(0.0f), fromPolar(alpha, beta));

	static float x0 = 0.5f;
	static float y0 = 0.8f;
	static float z0 = -2.0f;
	float x1 = x0 + exLen;
	float y1 = y0 + eyLen;

	glm::vec3 v[2][2];
	v[0][0] = glm::vec3(x0, y0, z0);
	v[0][1] = glm::vec3(x0, y1, z0);
	v[1][0] = glm::vec3(x1, y0, z0);
	v[1][1] = glm::vec3(x1, y1, z0);

	glm::vec3 n[4];
	auto calculate_n = [](glm::vec3 a, glm::vec3 b) {
		glm::vec3 c = glm::cross(a, b);
		return glm::normalize(c);
	};

	n[0] = calculate_n(v[0][0], v[1][0]);
	n[1] = calculate_n(v[1][0], v[1][1]);
	n[2] = calculate_n(v[1][1], v[0][1]);
	n[3] = calculate_n(v[0][1], v[0][0]);

	float rad = std::acos(glm::dot(-n[0], n[1]));

	ofSetColor(64);
	ofNoFill();
	ofDrawSphere(1.0);
	ofFill();
	
	for (int yi = 0; yi < 2; ++yi) {
		for (int xi = 0; xi < 2; ++xi) {
			glm::vec3 v_this = v[xi][yi];
			ofSetColor(ofColor::orange);
			ofDrawSphere(v_this, 0.01f);

			char text[128];
			sprintf(text, "v[%d, %d]", xi, yi);
			ofDrawBitmapString(text, v_this);
		}
	}

	{
		auto draw_projected = [](glm::vec3 a, glm::vec3 b) {
			int N = 30;
			ofPolyline line;
			for (int i = 0; i < N; ++i) {
				float s = (float)i / (N - 1);
				glm::vec3 p = glm::mix(a, b, s);
				line.addVertex(glm::normalize(p));
			}
			line.draw();

			ofDrawLine(a, glm::vec3(0.0f));
		};
		ofSetColor(ofColor::white);
		draw_projected(v[0][0], v[1][0]);
		draw_projected(v[1][0], v[1][1]);
		draw_projected(v[1][1], v[0][1]);
		draw_projected(v[0][1], v[0][0]);


		// normal
		auto draw_n = [](int i, glm::vec3 n, glm::vec3 va, glm::vec3 vb) {
			glm::vec3 q = glm::normalize(glm::normalize(va) + glm::normalize(vb));
			ofSetColor(ofColor::red);
			ofDrawLine(q, q + n * 0.2f);

			char text[128];
			sprintf(text, "n[%d]", i);
			ofDrawBitmapString(text, q);
		};
		draw_n(0, n[0], v[0][0], v[1][0]);
		draw_n(1, n[1], v[1][0], v[1][1]);
		draw_n(2, n[2], v[1][1], v[0][1]);
		draw_n(3, n[3], v[0][1], v[0][0]);
	}

	QuadPlane qp;
	qp.s = glm::dvec3(x0, y0, z0);
	qp.ex = glm::dvec3(exLen, 0.0, 0.0);
	qp.ey = glm::dvec3(0.0, eyLen, 0.0);
	SphereAreaSampler sampler(qp, glm::dvec3(0.0));

	{
		static rt::XoroshiroPlus128 random;
		ofMesh mesh;
		mesh.setMode(OF_PRIMITIVE_POINTS);
		for (int i = 0; i < 3000; ++i) {
			double u = random.uniform();
			double v = random.uniform();

			double phi_u = 0.0;
			// glm::dvec3 s = sampler.sample(u, v, &phi_u);
			glm::dvec3 s = sampler.sample(u, v);
			if (std::isfinite(s.x) == false ||
				std::isfinite(s.y) == false ||
				std::isfinite(s.z) == false) {
				printf("bad sample\n");
			}

			if (bang) {
				printf("[%d] %.10f\n",i, sin(phi_u));
			}

			mesh.addVertex(s);

			// projected
			mesh.addVertex(glm::normalize(s));
		}
		ofSetColor(255);
		mesh.draw();
	}
	bang = false;


	drawPlane(g_bugPlane);

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
	ImGui::Text("%.2f fps", ofGetFrameRate());
	ImGui::SliderFloat("x0", &x0, -5.0f, 5.0f);
	ImGui::SliderFloat("y0", &y0, -5.0f, 5.0f);
	ImGui::SliderFloat("z0", &z0, -5.0f, 5.0f);
	
	ImGui::SliderFloat("exLen", &exLen, 0.0f, 5.0f);
	ImGui::SliderFloat("eyLen", &eyLen, 0.0f, 5.0f);
	
	ImGui::SliderFloat("alpha (longitude)", &alpha, -glm::pi<float>(), glm::pi<float>());
	ImGui::SliderFloat("beta (polar)", &beta, 0.0f, glm::pi<float>());

	ImGui::Checkbox("bang", &bang);
	ImGui::Text("min phi(u) %.5f", sampler.minPhiU());
	ImGui::Text("max phi(u) %.5f", sampler.maxPhiU());

	ImGui::End();
}

//--------------------------------------------------------------
void ofApp::keyPressed(int key){
	if (key == ' ') {
		ofToggleFullscreen();
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
