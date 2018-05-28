#include "ofApp.h"
#include "bicubic.hpp"
#include "peseudo_random.hpp"
#include "microfacet.hpp"
#include "simpson.hpp"

#include <algorithm>
#include <tbb/tbb.h>

// zが上の座標系に移動する行列
inline glm::mat3 to_bxdf_basis_transform(const glm::vec3 &n) {
	glm::vec3 xaxis;
	glm::vec3 zaxis = n;
	glm::vec3 yaxis;
	if (0.999f < glm::abs(zaxis.z)) {
		xaxis = glm::normalize(glm::cross(glm::vec3(0.0f, 1.0f, 0.0f), zaxis));
	}
	else {
		xaxis = glm::normalize(glm::cross(glm::vec3(0.0f, 0.0f, 1.0f), zaxis));
	}
	yaxis = glm::cross(zaxis, xaxis);
	return glm::transpose(glm::mat3(xaxis, yaxis, zaxis));
}
inline glm::vec3 from_bxdf(const glm::vec3 &n, const glm::vec3 &bxdf_dir) {
	return glm::transpose(to_bxdf_basis_transform(n)) * bxdf_dir;
}

inline float evaluate_albedo(float o_theta, float alpha) {
	using namespace rt;
	glm::vec3 wo = rt::polar_to_cartesian(o_theta, 0.0f);

	int SampleCount = 100000;
	Xor *random = new Xor();

	glm::vec3 Ng(0.0f, 0.0f, 1.0f);

	double albedo_sum = 0.0;
	for (int i = 0; i < SampleCount; ++i) {
		float theta = std::atan(std::sqrt(-alpha * alpha * std::log(1.0f - random->uniform())));
		float phi = random->uniform(0.0f, glm::two_pi<double>());
		glm::vec3 sample = polar_to_cartesian(theta, phi);
		glm::vec3 harf = from_bxdf(Ng, sample);
		glm::vec3 wi = glm::reflect(-wo, harf);
		float pdf_omega = D_Beckman(Ng, harf, alpha) * glm::dot(Ng, harf) / (4.0f * glm::dot(wi, harf));
		if (glm::dot(Ng, wi) <= 0.0f) {
			continue;
		}

		glm::vec3 h = glm::normalize(wi + wo);
		float d = D_Beckman(Ng, h, alpha);
		float g = G2_height_correlated_beckmann(wi, wo, h, Ng, alpha);

		float cos_term_wo = glm::dot(Ng, wo);
		float cos_term_wi = glm::dot(Ng, wi);

		float brdf_without_f = d * g / (4.0f * cos_term_wo * cos_term_wi);

		float cosThetaFresnel = glm::dot(h, wo);
		float f = fresnel_dielectrics(cosThetaFresnel);

		float brdf = f * brdf_without_f * cos_term_wi;
		float albedo = brdf * cos_term_wi / pdf_omega;

		if (std::isfinite(albedo)) {
			albedo_sum += albedo;
		}

		//if (i % 500 == 0) {
		//	printf("%.8f\n", albedo_sum / (i + 1));
		//}
	}
	return albedo_sum / SampleCount;
}

class SpecularAlbedo {
public:
	void load(std::string path) {
		_specular_albedo.load(path);
	}
	float sample(float o_theta, float alpha) {
		float u = ofMap(alpha, 0.0001f, 1.0f, 0.0f, 1.0f);
		float v = ofMap(o_theta, 0.0001f, glm::radians(90.0f), 0.0f, 1.0f);
		float *p = _specular_albedo.getPixels().getPixels();
		float value = bicubic_2d(u, v, _specular_albedo.getWidth(), _specular_albedo.getHeight(), [&](int x, int y) { return p[y * (int)_specular_albedo.getWidth() + x]; });
		return value;
	}
	ofFloatImage _specular_albedo;
};
class SpecularAlbedoAvg {
public:
	void load(std::string path) {
		_specular_albedo_avg.load(path);
	}
	float sample(float alpha) {
		float u = ofMap(alpha, 0.0001f, 1.0f, 0.0f, 1.0f);
		float *p = _specular_albedo_avg.getPixels().getPixels();
		float value = bicubic_1d(u, _specular_albedo_avg.getWidth(), [&](int x) { return p[x]; });
		return value;
	}
	ofFloatImage _specular_albedo_avg;
};

//--------------------------------------------------------------
void ofApp::setup(){
	_camera.setNearClip(0.1f);
	_camera.setFarClip(100.0f);
	_camera.setDistance(5.0f);

	using namespace rt;

	//ofFloatImage image;
	//image.allocate(128, 128, OF_IMAGE_GRAYSCALE);
	//float *p = image.getPixels().getPixels();
	//for (int y = 0; y < 128; ++y) {
	//	for (int x = 0; x < 128; ++x) {
	//		p[y * 128 + x] = (float)x / 127;
	//	}
	//}
	// using namespace rt;

	// evaluate_albedo(glm::radians(30.0f), 0.2f);

	//int ComputeSize = 128;
	//ofFloatImage image;
	//image.allocate(ComputeSize, ComputeSize, OF_IMAGE_GRAYSCALE);
	//float *p = image.getPixels().getPixels();

	//tbb::parallel_for(tbb::blocked_range<int>(0, ComputeSize), [&](const tbb::blocked_range<int> &range) {
	//	for (int y = range.begin(); y < range.end(); ++y) {
	//		for (int x = 0; x < ComputeSize; ++x) {
	//			float alpha = ofMap(x, 0, ComputeSize - 1, 0.0001f, 1.0f);
	//			float o_theta = ofMap(y, 0, ComputeSize - 1, 0.0001f, glm::radians(90.0f));
	//			float albedo = evaluate_albedo(o_theta, alpha);
	//			p[y * ComputeSize + x] = albedo;
	//		}

	//		printf("%d line done.\n", y);
	//	}
	//});
	//image.save("specular_albedo.exr");

	SpecularAlbedo specularAlbedo;
	specularAlbedo.load("specular_albedo.exr");

	int Height = 10;
	int ComputeSize = 128;
	ofFloatImage image;
	image.allocate(ComputeSize, Height, OF_IMAGE_GRAYSCALE);
	float *p = image.getPixels().getPixels();
	for (int x = 0; x < ComputeSize; ++x) {
		float alpha = ofMap(x, 0, ComputeSize - 1, 0.0001f, 1.0f);

		double value = integrate_composite_simpson<100>([&](double phi) {
			return integrate_composite_simpson<100>([&](double theta) {
				double jacobian = std::sin(theta);
				glm::vec3 sample_m = polar_to_cartesian(theta, phi);
				double cosTheta = sample_m.z;
				double value = specularAlbedo.sample(theta, alpha) * cosTheta;
				return value * jacobian;
			}, 0.0, glm::pi<double>() * 0.5);
		}, 0.0, 2.0 * glm::pi<double>());

		double avg = (1.0 / glm::pi<double>()) * value;
		
		for (int y = 0; y < Height; ++y) {
			p[y * ComputeSize + x] = avg;
		}
	}

	image.save("albedo_avg.exr");

	//for (int y = 0; y < ComputeSize; ++y) {
	//	for (int x = 0; x < ComputeSize; ++x) {
	//		float alpha   = ofMap(x, 0, ComputeSize - 1, 0.0001f, 1.0f);
	//		float o_theta = ofMap(x, 0, ComputeSize - 1, 0.0001f, glm::radians(90.0f));
	//		float albedo = evaluate_albedo(o_theta, alpha);
	//		p[y * 128 + x] = albedo;
	//	}
	//	printf("%d line done.\n", y);
	//}
	

	//float alpha = 0.2f;
	//glm::vec3 wo = rt::polar_to_cartesian(glm::radians(30.0f), 0.0f);
	//int SampleCount = 100000;
	//Xor *random = new Xor();

	//glm::vec3 Ng(0.0f, 0.0f, 1.0f);

	//float albedo_sum = 0.0f;
	//for (int i = 0; i < SampleCount; ++i) {
	//	float theta = std::atan(std::sqrt(-alpha * alpha * std::log(1.0f - random->uniform())));
	//	float phi = random->uniform(0.0f, glm::two_pi<double>());
	//	glm::vec3 sample = polar_to_cartesian(theta, phi);
	//	glm::vec3 harf = from_bxdf(Ng, sample);
	//	glm::vec3 wi = glm::reflect(-wo, harf);
	//	float pdf_omega = D_Beckman(Ng, harf, alpha) * glm::dot(Ng, harf) / (4.0f * glm::dot(wi, harf));
	//	if (glm::dot(Ng, wi) <= 0.0f) {
	//		continue;
	//	}

	//	glm::vec3 h = glm::normalize(wi + wo);
	//	float d = D_Beckman(Ng, h, alpha);
	//	float g = G2_height_correlated_beckmann(wi, wo, h, Ng, alpha);

	//	float cos_term_wo = glm::dot(Ng, wo);
	//	float cos_term_wi = glm::dot(Ng, wi);

	//	float brdf_without_f = d * g / (4.0f * cos_term_wo * cos_term_wi);

	//	float cosThetaFresnel = glm::dot(h, wo);
	//	float f = fresnel_dielectrics(cosThetaFresnel);

	//	float brdf = f * brdf_without_f * cos_term_wi;

	//	albedo_sum += brdf * cos_term_wi / pdf_omega;

	//	if (i % 500 == 0) {
	//		printf("%.8f\n", albedo_sum / (i + 1));
	//	}
	//}
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

	int N = 50;
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
