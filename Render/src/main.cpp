#include "alembic_loader.hpp"
#include "integrator.hpp"
#include "online.hpp"

#include <functional>
#include <random>
#include <xmmintrin.h>
#include <pmmintrin.h>

#include "ofMain.h"

static const double kRenderTime = 123.0;
// static const double kRenderTime = 60 * 60;

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

// step(cur frame), save(spp)
inline void render(rt::Stopwatch *main_sw, std::function<void(int)> step, std::function<void(int)> save, double duration, double save_interval) {
	double safety_duration = 1.0;

	rt::Stopwatch save_interval_time;

	rt::OnlineMean<double> step_duration;
	for (int i = 0; ; ++i) {
		rt::Stopwatch sw;
		step(i);
		step_duration.addSample(sw.elapsed());

		if (duration - safety_duration < main_sw->elapsed() + step_duration.mean()) {
			save(i + 1);
			break;
		}

		if (save_interval < save_interval_time.elapsed()) {
			save(i + 1);
			save_interval_time = rt::Stopwatch();
		}
	}
}


//========================================================================
int main( ){
	rt::Stopwatch sw;

	rt::CoupledBRDFConductor::load(
		ofToDataPath("baked/albedo_specular_conductor.bin").c_str(),
		ofToDataPath("baked/albedo_specular_conductor_avg.bin").c_str());
	rt::CoupledBRDFDielectrics::load(
		ofToDataPath("baked/albedo_specular_dielectrics.bin").c_str(),
		ofToDataPath("baked/albedo_specular_dielectrics_avg.bin").c_str());

	rt::CoupledBRDFVelvet::load(
		ofToDataPath("baked/albedo_velvet.bin").c_str(),
		ofToDataPath("baked/albedo_velvet_avg.bin").c_str());
	
	std::shared_ptr<rt::Scene> scene(new rt::Scene());
	rt::loadFromABC(ofToDataPath("cornelbox.abc").c_str(), *scene);
	printf("setup %f seconds\n", sw.elapsed());

	std::shared_ptr<rt::PTRenderer> renderer(new rt::PTRenderer(scene));

	render(&sw, [&](int frame) {
		renderer->step();
		printf("frame %03d, %.1f sec\n", frame, sw.elapsed());
	}, [&](int spp) {
		ofPixels image = toOf(renderer->_image);

		char name[128];
		sprintf(name, "../../rendered_images/image_%d_spp.png", spp);
		ofSaveImage(image, name);
		printf("save as %s\n", name);
	}, kRenderTime, 15.0);
}
