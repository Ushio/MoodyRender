#include "ofApp.h"
#include "bicubic.hpp"
#include "peseudo_random.hpp"
#include "microfacet.hpp"
#include "simpson.hpp"

#include <algorithm>
#include <tbb/tbb.h>

// zが上の座標系に移動する行列
inline glm::mat3 to_bxdf_basis_transform(const glm::dvec3 &n) {
	glm::dvec3 xaxis;
	glm::dvec3 zaxis = n;
	glm::dvec3 yaxis;
	if (0.999 < glm::abs(zaxis.z)) {
		xaxis = glm::normalize(glm::cross(glm::dvec3(0.0, 1.0, 0.0), zaxis));
	}
	else {
		xaxis = glm::normalize(glm::cross(glm::dvec3(0.0, 0.0, 1.0), zaxis));
	}
	yaxis = glm::cross(zaxis, xaxis);
	return glm::transpose(glm::mat3(xaxis, yaxis, zaxis));
}
inline glm::dvec3 from_bxdf(const glm::dvec3 &n, const glm::dvec3 &bxdf_dir) {
	return glm::transpose(to_bxdf_basis_transform(n)) * bxdf_dir;
}

inline double evaluate_albedo(double o_theta, double alpha) {
	using namespace rt;
	glm::dvec3 wo = rt::polar_to_cartesian(o_theta, 0.0);

	int SampleCount = 100000;
	Xor *random = new Xor();

	glm::dvec3 Ng(0.0, 0.0, 1.0);

	double albedo_sum = 0.0;
	for (int i = 0; i < SampleCount; ++i) {
		double theta = std::atan(std::sqrt(-alpha * alpha * std::log(1.0 - random->uniform())));
		double phi = random->uniform(0.0, glm::two_pi<double>());
		glm::dvec3 sample = polar_to_cartesian(theta, phi);
		glm::dvec3 harf = from_bxdf(Ng, sample);
		glm::dvec3 wi = glm::reflect(-wo, harf);
		double pdf_omega = D_Beckman(Ng, harf, alpha) * glm::dot(Ng, harf) / (4.0 * glm::dot(wi, harf));
		if (glm::dot(Ng, wi) <= 0.0) {
			continue;
		}

		glm::dvec3 h = glm::normalize(wi + wo);
		double d = D_Beckman(Ng, h, alpha);
		double g = G2_height_correlated_beckmann(wi, wo, h, Ng, alpha);

		double cos_term_wo = glm::dot(Ng, wo);
		double cos_term_wi = glm::dot(Ng, wi);

		double brdf_without_f = d * g / (4.0 * cos_term_wo * cos_term_wi);

		double cosThetaFresnel = glm::dot(h, wo);
		double f = fresnel_dielectrics(cosThetaFresnel);

		double brdf = f * brdf_without_f * cos_term_wi;
		double albedo = brdf * cos_term_wi / pdf_omega;

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
	double sample(double o_theta, double alpha) {
		double u = ofMap(alpha, 0.0001, 1.0, 0.0, 1.0);
		double v = ofMap(o_theta, 0.0001, glm::radians(90.0), 0.0, 1.0);
		double *p = _specular_albedo.getPixels().getPixels();
		double value = bicubic_2d(u, v, _specular_albedo.getWidth(), _specular_albedo.getHeight(), [&](int x, int y) { return p[y * (int)_specular_albedo.getWidth() + x]; });
		return value;
	}
	ofFloatImage _specular_albedo;
};
class SpecularAlbedoAvg {
public:
	void load(std::string path) {
		_specular_albedo_avg.load(path);
	}
	double sample(double alpha) {
		double u = ofMap(alpha, 0.0001, 1.0, 0.0, 1.0);
		double *p = _specular_albedo_avg.getPixels().getPixels();
		double value = bicubic_1d(u, _specular_albedo_avg.getWidth(), [&](int x) { return p[x]; });
		return value;
	}
	ofFloatImage _specular_albedo_avg;
};

//--------------------------------------------------------------
void ofApp::setup(){
	_camera.setNearClip(0.1);
	_camera.setFarClip(100.0);
	_camera.setDistance(5.0);

	using namespace rt;

	//ofFloatImage image;
	//image.allocate(128, 128, OF_IMAGE_GRAYSCALE);
	//double *p = image.getPixels().getPixels();
	//for (int y = 0; y < 128; ++y) {
	//	for (int x = 0; x < 128; ++x) {
	//		p[y * 128 + x] = (double)x / 127;
	//	}
	//}
	// using namespace rt;

	// evaluate_albedo(glm::radians(30.0), 0.2);

	//int ComputeSize = 128;
	//ofFloatImage image;
	//image.allocate(ComputeSize, ComputeSize, OF_IMAGE_GRAYSCALE);
	//double *p = image.getPixels().getPixels();

	//tbb::parallel_for(tbb::blocked_range<int>(0, ComputeSize), [&](const tbb::blocked_range<int> &range) {
	//	for (int y = range.begin(); y < range.end(); ++y) {
	//		for (int x = 0; x < ComputeSize; ++x) {
	//			double alpha = ofMap(x, 0, ComputeSize - 1, 0.0001, 1.0);
	//			double o_theta = ofMap(y, 0, ComputeSize - 1, 0.0001, glm::radians(90.0));
	//			double albedo = evaluate_albedo(o_theta, alpha);
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
	double *p = image.getPixels().getPixels();
	for (int x = 0; x < ComputeSize; ++x) {
		double alpha = ofMap(x, 0, ComputeSize - 1, 0.0001, 1.0);

		double value = integrate_composite_simpson<100>([&](double phi) {
			return integrate_composite_simpson<100>([&](double theta) {
				double jacobian = std::sin(theta);
				glm::dvec3 sample_m = polar_to_cartesian(theta, phi);
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
	//		double alpha   = ofMap(x, 0, ComputeSize - 1, 0.0001, 1.0);
	//		double o_theta = ofMap(x, 0, ComputeSize - 1, 0.0001, glm::radians(90.0));
	//		double albedo = evaluate_albedo(o_theta, alpha);
	//		p[y * 128 + x] = albedo;
	//	}
	//	printf("%d line done.\n", y);
	//}
	

	//double alpha = 0.2;
	//glm::dvec3 wo = rt::polar_to_cartesian(glm::radians(30.0), 0.0);
	//int SampleCount = 100000;
	//Xor *random = new Xor();

	//glm::dvec3 Ng(0.0, 0.0, 1.0);

	//double albedo_sum = 0.0;
	//for (int i = 0; i < SampleCount; ++i) {
	//	double theta = std::atan(std::sqrt(-alpha * alpha * std::log(1.0 - random->uniform())));
	//	double phi = random->uniform(0.0, glm::two_pi<double>());
	//	glm::dvec3 sample = polar_to_cartesian(theta, phi);
	//	glm::dvec3 harf = from_bxdf(Ng, sample);
	//	glm::dvec3 wi = glm::reflect(-wo, harf);
	//	double pdf_omega = D_Beckman(Ng, harf, alpha) * glm::dot(Ng, harf) / (4.0 * glm::dot(wi, harf));
	//	if (glm::dot(Ng, wi) <= 0.0) {
	//		continue;
	//	}

	//	glm::dvec3 h = glm::normalize(wi + wo);
	//	double d = D_Beckman(Ng, h, alpha);
	//	double g = G2_height_correlated_beckmann(wi, wo, h, Ng, alpha);

	//	double cos_term_wo = glm::dot(Ng, wo);
	//	double cos_term_wi = glm::dot(Ng, wi);

	//	double brdf_without_f = d * g / (4.0 * cos_term_wo * cos_term_wi);

	//	double cosThetaFresnel = glm::dot(h, wo);
	//	double f = fresnel_dielectrics(cosThetaFresnel);

	//	double brdf = f * brdf_without_f * cos_term_wi;

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
	ofRotateZ(90.0);
	ofSetColor(128);
	ofDrawGridPlane(1.0);
	ofPopMatrix();

	ofPushMatrix();
	ofDrawAxis(50);
	ofPopMatrix();

#if 0
	rt::Xor random;
	double image[4];
	for (int i = 0; i < 4; ++i) {
		image[i] = random.uniform();
	}

	ofSetColor(255);
	for (int i = 0; i < 4; ++i) {
		ofDrawSphere(i, image[i], 0.05);
	}

	ofMesh mesh;
	mesh.setMode(OF_PRIMITIVE_POINTS);

	int N = 1000; 
	for (int i = 0; i < N; ++i) {
		double x = ofMap(i, 0, N - 1, -0.2, 1.2);
		double value = bicubic_1d(x, 4, [&](int x) { return image[x]; });
		mesh.addVertex(glm::dvec3(x * 3, value, 0.0));
	}
	mesh.draw();
#else
	const int X_SIZE = 6;
	const int Y_SIZE = 6;

	rt::Xor random;
	double image[Y_SIZE][X_SIZE];
	for (int z = 0; z < Y_SIZE; ++z) {
		for (int x = 0; x < X_SIZE; ++x) {
			image[z][x] = random.uniform();
		}
	}

	ofSetColor(255);
	for (int z = 0; z < Y_SIZE; ++z) {
		for (int x = 0; x < X_SIZE; ++x) {
			ofDrawSphere(x, image[z][x], z, 0.05);
		}
	}

	ofMesh mesh;
	mesh.setMode(OF_PRIMITIVE_POINTS);

	int N = 50;
	for (int yi = 0; yi < N; ++yi) {
		double y = ofMap(yi, 0, N - 1, -0.5, 1.5);

		for (int xi = 0; xi < N; ++xi) {
			double x = ofMap(xi, 0, N - 1, -0.5, 1.5);

			double value = bicubic_2d(x, y, X_SIZE, Y_SIZE, [&](int x, int y) { return image[y][x]; });
			mesh.addVertex(glm::dvec3(x * (X_SIZE - 1), value, y * (Y_SIZE - 1)));
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
