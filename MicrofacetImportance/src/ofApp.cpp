#include "ofApp.h"
#include "ofxImGuiLite.hpp"

#include "peseudo_random.hpp"
#include "microfacet.hpp"
#include "material.hpp"
#include "simpson_helper.hpp"

#include "randomsampler.hpp"

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

using namespace rt;
inline glm::dvec3 sample_on_unit_sphere_super(PeseudoRandom *random) {
	double x1;
	double x2;
	double S;
	do {
		x1 = random->uniform(-1.0, 1.0);
		x2 = random->uniform(-1.0, 1.0);
		S = x1 * x1 + x2 * x2;
	} while (S >= 1.0);

	double two_sqrt_one_minus_s = 2.0 * std::sqrt(std::max(1.0 - S, 0.0));
	return glm::dvec3(
		x1 * two_sqrt_one_minus_s,
		x2 * two_sqrt_one_minus_s,
		1.0 - 2.0 * S);
}

//glm::dvec3 sample_on_unit_sphere_super(PeseudoRandom *random) {
//	double x1;
//	double x2;
//	double S;
//	do {
//		// x1 = random->uniform(-1.0, 1.0);
//		// x2 = random->uniform(-1.0, 1.0);
//		static const double kSqrt2 = std::sqrt(2.0);
//		static const double kSqrt3 = std::sqrt(3.0);
//
//		x1 = random->uniform(-1.0, 1.0) * kSqrt3;
//		x2 = random->uniform(-1.0, 1.0);
//		
//		//double D = -kSqrt3 * std::fabs(x1) + 2.0;
//		//if (D < std::fabs(x2)) {
//		//	x1 = (x1 - std::copysign(kSqrt3, x1));
//		//	x2 = (x2 - std::copysign(1.0, x2));
//		//}
//
//		S = x1 * x1 + x2 * x2;
//	} while (1.0 <= S);
//
//	double two_sqrt_one_minus_s = 2.0 * std::sqrt(std::max(1.0 - S, 0.0));
//	return glm::dvec3(
//		x1 * two_sqrt_one_minus_s,
//		x2 * two_sqrt_one_minus_s,
//		1.0 - 2.0 * S);
//}

//glm::dvec3 sample_on_unit_sphere_super(PeseudoRandom *random) {
//	double x1;
//	double x2;
//	double S;
//	do {
//		x1 = random->uniform(-1.0, 1.0);
//		x2 = random->uniform(-1.0, 1.0);
//
//		S = x1 * x1 + x2 * x2;
//
//		if (S < 1.0) {
//			break;
//		}
//
//		//static const double r = 3.0f - 2.0 * sqrt(2.0);
//		//static const double r_inv = 1.0 / r;
//		//static const double c_xy = 1.0f - r;
//		//static const double K = c_xy * r_inv;
//
//		//x1 = std::abs(x1);
//		//x2 = std::abs(x2);
//
//		// x1 = (x1 - c_xy) * r_inv;
//		// x2 = (x2 - c_xy) * r_inv;
//		// x1 = x1 * r_inv - K;
//		//x1 = std::fma(std::abs(x1), r_inv, -K);
//		//x2 = std::fma(std::abs(x2), r_inv, -K);
//		//S = x1 * x1 + x2 * x2;
//
//		static const double r = sqrt(2.0) - 1.0;
//		static const double r_inv = 1.0 / r;
//		//x1 = (x1 - std::copysign(1.0, x1)) * r_inv;
//		//x2 = (x2 - std::copysign(1.0, x2)) * r_inv;
//		x1 = x1 * r_inv - std::copysign(r_inv, x1);
//		x2 = x2 * r_inv - std::copysign(r_inv, x2);
//		S = x1 * x1 + x2 * x2;
//		if (S < 1.0) {
//			break;
//		}
//	} while (true);
//
//	double two_sqrt_one_minus_s = 2.0 * std::sqrt(std::max(1.0 - S, 0.0));
//	return glm::dvec3(
//		x1 * two_sqrt_one_minus_s,
//		x2 * two_sqrt_one_minus_s,
//		1.0 - 2.0 * S);
//}

// Marsaglia_64f
//inline glm::dvec3 sample_on_unit_sphere(PeseudoRandom *random, std::vector<double> *values) {
//	double x1;
//	double x2;
//	double S;
//	do {
//		x1 = random->uniform(-1.0, 1.0);
//		x2 = random->uniform(-1.0, 1.0);
//		values->push_back((x1 + 1.0) * 0.5);
//		values->push_back((x2 + 1.0) * 0.5);
//
//		S = x1 * x1 + x2 * x2;
//	} while (S >= 1.0);
//
//	double two_sqrt_one_minus_s = 2.0 * std::sqrt(std::max(1.0 - S, 0.0));
//	return glm::dvec3(
//		x1 * two_sqrt_one_minus_s,
//		x2 * two_sqrt_one_minus_s,
//		1.0 - 2.0 * S);
//}

// PBRT_64f.txt
//inline glm::dvec3 sample_on_unit_sphere(PeseudoRandom *random, std::vector<double> *values) {
//	double z = random->uniform(-1.0, 1.0);
//	double phi = random->uniform(0.0, glm::two_pi<double>());
//	values->push_back((z + 1.0) * 0.5);
//	values->push_back(phi / glm::two_pi<double>());
//
//	double sq_one_minus_zz = std::sqrt(1.0 - z * z);
//	double x = sq_one_minus_zz * std::cos(phi);
//	double y = sq_one_minus_zz * std::sin(phi);
//	return glm::dvec3(x, y, z);
//}

//inline glm::dvec3 sample_on_unit_sphere(rt::PeseudoRandom *random) {
//	glm::dvec3 d;
//	double sq = 0.0;
//	do {
//		d.x = random->uniform(-1.0, 1.0);
//		d.y = random->uniform(-1.0, 1.0);
//		d.z = random->uniform(-1.0, 1.0);
//
//		sq = glm::length2(d);
//	} while (sq < 1.0e-9 || 1.0 < sq);
//	d /= glm::sqrt(sq);
//	return d;
//}

//inline glm::vec3 sample_on_unit_sphere(PeseudoRandom *random) {
//	float x1;
//	float x2;
//	float S;
//	do {
//		x1 = random->uniformf(-1.0f, 1.0f);
//		x2 = random->uniformf(-1.0f, 1.0f);
//		S = x1 * x1 + x2 * x2;
//	} while (S >= 1.0f);
//
//	float two_sqrt_one_minus_s = 2.0f * sqrt(1.0f - S);
//	return glm::vec3(
//		x1 * two_sqrt_one_minus_s,
//		x2 * two_sqrt_one_minus_s,
//		1.0f - 2.0f * S);
//}

//inline glm::vec3 sample_on_unit_sphere(PeseudoRandom *random) {
//	float z = random->uniformf(-1.0f, 1.0f);
//	float phi = random->uniformf(0.0f, glm::two_pi<float>());
//	float sq_one_minus_zz = std::sqrt(1.0f - z * z);
//	float x = sq_one_minus_zz * std::cos(phi);
//	float y = sq_one_minus_zz * std::sin(phi);
//	return glm::vec3(x, y, z);
//}

template <class RandomType>
void plot(const char *filename, std::function<double(PeseudoRandom *)> z) {
	FILE *fp = fopen(filename, "w");
	printf("%s\n", filename);
	fprintf(fp, "%s\n", filename);

	RandomType random;
	
	for (int i = 0; i < 10000; ++i) {
		fprintf(fp, "%.15f\n", z(&random));
	}

	fclose(fp);
}


template <class Real>
class OnlineMean {
public:
	void addSample(Real x) {
		_sum += x;
		_n++;
	}
	Real mean() const {
		return _sum / _n;
	}
	int sampleCount() const {
		return _n;
	}
private:
	int _n = 0;
	Real _sum = Real(0.0);
};

class SphereStore {
public:
	SphereStore(int w, int h):_w(w), _h(h), _storage(w * h){

	}
	void addSample(double theta, double phi, double value) {
		int x = xFromTheta(theta);
		int y = yFromPhi(phi);
		_storage[y * _w + x].addSample(value);
	}
	ofFloatImage toImage() const {
		ofFloatImage image;
		image.allocate(_w, _h, OF_IMAGE_GRAYSCALE);
		float *p = image.getPixels().getPixels();
		for (int i = 0; i < _storage.size(); ++i) {
			p[i] = 0 < _storage[i].sampleCount() ? _storage[i].mean() : 0;
		}
		return image;
	}
private:
	// 0 ~ pi
	int xFromTheta(double theta) const {
		double xf = theta / glm::pi<double>();
		int x = (int)std::floor(xf * _w);
		x = glm::clamp(x, 0, _w - 1);
		return x;
	}
	int yFromPhi(double phi) const {
		double yf = phi / glm::two_pi<double>();
		int y = (int)std::floor(yf * _h);
		y = glm::clamp(y, 0, _h - 1);
		return y;
	}
	int _w = 0;
	int _h = 0;
	std::vector<OnlineMean<double>> _storage;
};

//--------------------------------------------------------------
void ofApp::setup(){
	ofxImGuiLite::initialize();

	using namespace rt;
	/*
	XoroshiroPlus128 random;
	std::vector<double> values;

	SphereStore store(128, 128);
	for (int i = 0; i < 10000000; ++i) {
		values.clear();
		glm::dvec3 sample = sample_on_unit_sphere(&random, &values);
		double theta = std::acos(sample.z);
		double phi = glm::pi<double>() + std::copysign(1.0, sample.y) * std::acos(sample.x / std::sqrt(sample.x * sample.x + sample.y * sample.y));
	
		for (double value : values) {
			store.addSample(theta, phi, value);
		}
		// store.addSample(theta, phi, theta / glm::pi<double>());
		// store.addSample(theta, phi, phi / glm::two_pi<double>());
		// store.addSample(theta, phi, 0.5);
	}
	ofFloatImage image = store.toImage();
	// image.save("image.exr");
	ofImage(image).save("image.png");
	*/

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

	ofSetCircleResolution(100);

	ofSetColor(255);
	ofNoFill();
	ofDrawCircle(0, 0, 1);
	ofDrawCircle(sqrt(3.0), 1, 1);
	ofDrawCircle(sqrt(3.0), -1, 1);
	ofFill();

	// float xc = (2 - sqrt(2.0f));
	// ofDrawCircle(1, 1 - (2 - sqrt(2.0f)), 0.01);
	// ofDrawCircle(1 - (2 - sqrt(2.0f)), 1, 0.01);
	//float r = 3.0f - 2.0 * sqrt(2.0);
	//ofDrawCircle(1 - r, 1 - r, r);
	//ofFill();

	//float x = ofMap(ofGetMouseX(), 0, ofGetWidth(), -1, 1);
	//float y = ofMap(ofGetMouseY(), 0, ofGetHeight(), 1, -1);

	//ofSetColor(255);
	//ofNoFill();
	//ofDrawCircle(0, 0, 1.0);
	//ofFill();

	// ofDrawLine(0.0, -0.0 + sqrt(2.0), 1.0, -1.0 + sqrt(2.0));

	//float S = x * x + y * y;
	//if (S < 1.0f) {
	//	ofDrawSphere(x, y, 0.0, 0.05);
	//}
	//else {
	//	ofDrawSphere(x - std::copysign(1.0, x), y - std::copysign(1.0, y), 0.0, 0.05);
	//}

	//{
	//	Xor64 random;
	//	ofMesh mesh;
	//	mesh.setMode(OF_PRIMITIVE_POINTS);
	//	for (int i = 0; i < 10000; ++i) {
	//		static const double kSqrt2 = std::sqrt(2.0);
	//		static const double kSqrt3 = std::sqrt(3.0);

	//		float x = random.uniformf(-1, 1) * kSqrt3;
	//		float y = random.uniformf(-1, 1);
	//		float S = x * x + y * y;
	//		if (S < 1.0f) {

	//		}
	//		else {
	//			//static const double r = sqrt(2.0) - 1.0;
	//			//static const double r_inv = 1.0 / r;
	//			//x = (x - std::copysign(1.0, x)) * r_inv;
	//			//y = (y - std::copysign(1.0, y)) * r_inv;

	//			//mesh.addVertex(glm::vec3(x, y, 0.0f));

	//			double D = -kSqrt3 * std::fabs(x) + 2.0;
	//			if (D < std::fabs(y)) {
	//				x = (x - std::copysign(kSqrt3, x));
	//				y = (y - std::copysign(1.0, y));
	//			}
	//			mesh.addVertex(glm::vec3(x, y, 0.0f));
	//		}
	//	}
	//	mesh.draw();
	//}

	//{
	//	Xor64 random;
	//	ofMesh mesh;
	//	mesh.setMode(OF_PRIMITIVE_POINTS);
	//	for (int i = 0; i < 10000; ++i) {
	//		float x = random.uniformf(-1, 1);
	//		float y = random.uniformf(-1, 1);
	//		float S = x * x + y * y;
	//		if (S < 1.0f) {

	//		}
	//		else {
	//			x = x - std::copysign(1.0f, x);
	//			y = y - std::copysign(1.0f, y);
	//			static const float T = 1.0f / std::sqrt(2.0f);
	//			static const float scale = 1.0f / (std::sqrt(2.0f) - 1.0f);
	//			static const float K = T * scale;
	//			//float rx = x * T - y * T;
	//			//float ry = x * T + y * T;
	//			//mesh.addVertex(glm::vec3(rx * scale, ry * scale, 0.0f));
	//			float rx = (x - y) * K;
	//			float ry = (x + y) * K;
	//			mesh.addVertex(glm::vec3(rx, ry, 0.0f));
	//		}
	//	}
	//	mesh.draw();
	//}

	// 

	//double cosTheta = cos(theta);
	//double sinTheta = sin(theta);
	//glm::dvec3 wo = glm::dvec3(sinTheta, 0.0, cosTheta);
	//glm::dvec3 Ng(0.0, 0.0, 1.0);

	//ofSetColor(255, 0, 0);
	//ofDrawLine(glm::vec3(), glm::vec3(wo.x, wo.z, wo.y));

	//{
	//	ofSetColor(255);

	//	auto geoffrey = [](float t) {
	//		glm::vec3 r = glm::vec3(t * 2.0f) - glm::vec3(1.8, 1.14f, 0.3f);
	//		return glm::clamp(glm::vec3(1.0f) - r * r, glm::vec3(0.0), glm::vec3(1.0));
	//	};
	//	ofMesh mesh;
	//	mesh.setMode(OF_PRIMITIVE_POINTS);
	//	XoroshiroPlus128 random;
	//	std::vector<double> values;
	//	for (int i = 0; i < 10000; ++i) {
	//		values.clear();
	//		glm::dvec3 sample = sample_on_unit_sphere(&random, &values);
	//		double theta = std::acos(sample.z);
	//		double phi = glm::pi<double>() + std::copysign(1.0, sample.y) * std::acos(sample.x / std::sqrt(sample.x * sample.x + sample.y * sample.y));

	//		mesh.addVertex(glm::dvec3(sample.x, sample.z, sample.y));

	//		// glm::vec3 c = geoffrey(phi / glm::two_pi<double>());
	//		glm::vec3 c = geoffrey(theta / glm::pi<double>());
	//		mesh.addColor(ofColor(255 * c.x, 255 * c.y, 255 * c.z));

	//		//char text[256];
	//		//// sprintf(text, "%.2f, %.2f", theta, phi);
	//		//sprintf(text, "%.2f", phi);
	//		//ofDrawBitmapString(text, sample.x, sample.z, sample.y);
	//	}
	//	
	//	mesh.draw();
	//}

// bench
	{
		XoroshiroPlus128 random;
		ofMesh mesh;
		mesh.setMode(OF_PRIMITIVE_POINTS);
		glm::dvec3 s = glm::dvec3(0.0);
		for (int i = 0; i < 3000000; ++i) {
			glm::dvec3 wi = sample_on_unit_sphere_super(&random);
			if (i % 100 == 0) {
				mesh.addVertex(glm::dvec3(wi.x, wi.z, wi.y));
			}
			s += wi;
		}
		ofSetColor(255);
		mesh.draw();
		printf("%f, %f, %f\n", s.x, s.y, s.z);
	}

	//{
	//	ofMesh mesh;
	//	mesh.setMode(OF_PRIMITIVE_POINTS);
	//	for (int i = 0; i < 3000; ++i) {
	//		glm::dvec3 wi = VelvetSampler::sample(&random, alpha, wo, Ng);
	//		mesh.addVertex(glm::dvec3(wi.x, wi.z, wi.y));
	//	}
	//	ofSetColor(255);
	//	mesh.draw();

	//}
	//{
	//	ofPolyline line;
	//	int N = 3000;
	//	for (int i = 0; i < N; ++i) {
	//		double theta = ofMap(i, 0, N - 1, 0, glm::two_pi<double>());
	//		double cosTheta = cos(theta);
	//		double sinTheta = sin(theta);
	//		glm::dvec3 wi = glm::dvec3(sinTheta, 0.0, cosTheta);
	//		// glm::dvec3 half = glm::normalize(wi + wo);
	//		double p = D_velvet(Ng, wi, alpha);

	//		line.addVertex(glm::dvec3(wi.x, wi.z, wi.y) * p);
	//	}
	//	line.draw();
	//}

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
	ImGui::Text("%.2f fps", ofGetFrameRate());
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
