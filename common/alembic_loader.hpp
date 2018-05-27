#pragma once

#include <Alembic/Abc/All.h>
#include <Alembic/AbcCoreOgawa/All.h>
#include <Alembic/AbcGeom/All.h>

#include "render_object.hpp"

#include <functional>

namespace rt {
	using namespace Alembic::Abc;
	using namespace Alembic::AbcGeom;

	inline void visitProperties(ICompoundProperty props,
		std::function<void(IScalarProperty, std::string)> scalarProp,
		std::function<void(IArrayProperty, std::string)> arrayProp,
		std::function<void(ICompoundProperty, std::string)> compoundProp,
		std::string dir = "") {

		for (int i = 0; i < props.getNumProperties(); ++i) {
			auto header = props.getPropertyHeader(i);
			auto child_header = props.getPropertyHeader(i);
			auto name = dir + "/" + child_header.getName();

			auto child_type = child_header.getPropertyType();
			if (child_type == AbcA::PropertyType::kScalarProperty) {
				auto child_prop = IScalarProperty(props, child_header.getName());
				scalarProp(child_prop, name);
			}
			if (child_type == AbcA::PropertyType::kArrayProperty) {
				auto child_prop = IArrayProperty(props, child_header.getName());
				arrayProp(child_prop, name);
			}
			if (child_type == AbcA::PropertyType::kCompoundProperty) {
				auto child_prop = ICompoundProperty(props, child_header.getName());
				compoundProp(child_prop, name);
				visitProperties(child_prop, scalarProp, arrayProp, compoundProp, name);
			}
		}
	}

	inline void printProperties(ICompoundProperty props, std::string dir = "") {
		visitProperties(props,
			[](IScalarProperty prop, std::string name) { printf("  %s [kScalarProperty]\n", name.c_str()); },
			[](IArrayProperty prop, std::string name) { printf("  %s [IArrayProperty]\n", name.c_str()); },
			[](ICompoundProperty prop, std::string name) { printf("  %s [ICompoundProperty]\n", name.c_str()); }
		);
	}
	inline std::vector<string> propertyMaterial(ICompoundProperty props) throw(std::exception) {
		const char *vals_name = "/.geom/.arbGeomParams/Material/.vals";
		const char *indices_name = "/.geom/.arbGeomParams/Material/.indices";
		std::vector<std::string> vals;
		std::vector<uint32_t> indices;

		bool foundVals = false;
		bool foundIndices = false;

		visitProperties(props,
			[](IScalarProperty prop, std::string name) {},
			[&](IArrayProperty prop, std::string name) {
			if (name == indices_name) {
				foundIndices = true;

				Abc::IUInt32ArrayProperty arrayProp(prop.getParent(), prop.getName());
				UInt32ArraySamplePtr sample;
				arrayProp.get(sample);

				indices.resize(sample->size());
				for (int j = 0; j < sample->size(); ++j) {
					indices[j] = sample->get()[j];
				}
			}
			else if (name == vals_name) {
				foundVals = true;

				Abc::IStringArrayProperty stringProp(prop.getParent(), prop.getName());
				StringArraySamplePtr sample;
				stringProp.get(sample);

				vals.resize(sample->size());
				for (int j = 0; j < sample->size(); ++j) {
					vals[j] = sample->get()[j];
				}
			}
		},
			[](ICompoundProperty prop, std::string name) {}
		);

		if (foundVals == false || foundIndices == false) {
			throw std::exception("key not found");
		}

		std::vector<std::string> expands(indices.size());
		for (int i = 0; i < indices.size(); ++i) {
			expands[i] = vals[indices[i]];
		}

		return expands;
	}

	inline std::vector<glm::vec3> toVec3Array(V3fArraySamplePtr sample) {
		std::vector<glm::vec3> values(sample->size());
		for (int j = 0; j < sample->size(); ++j) {
			auto value = sample->get()[j];
			values[j] = glm::vec3(value.x, value.y, value.z);
		}
		return values;
	}

	inline std::vector<glm::vec3> propertyArrayVec3(ICompoundProperty props, const char *dir) throw(std::exception) {
		bool found = false;
		std::vector<glm::vec3> values;
		visitProperties(props,
			[](IScalarProperty prop, std::string name) {},
			[&](IArrayProperty prop, std::string name) {
			if (name == dir) {
				Abc::IV3fArrayProperty arrayProp(prop.getParent(), prop.getName());
				V3fArraySamplePtr sample;
				arrayProp.get(sample);

				found = true;
				values.resize(sample->size());
				for (int j = 0; j < sample->size(); ++j) {
					auto value = sample->get()[j];
					values[j] = glm::vec3(value.x, value.y, value.z);
				}
			}
		},
			[](ICompoundProperty prop, std::string name) {}
		);
		if (found == false) {
			throw std::exception("key not found");
		}
		return values;
	}
	inline float propertyScalarFloat(ICompoundProperty props, const char *dir) {
		bool found = false;
		float value = 0.0f;
		visitProperties(props,
			[&](IScalarProperty prop, std::string name) {
			if (name == dir) {
				found = true;
				IFloatProperty floatProp(prop.getParent(), prop.getName());
				floatProp.get(value);
			}
		},
			[](IArrayProperty prop, std::string name) {},
			[](ICompoundProperty prop, std::string name) {}
		);
		if (found == false) {
			throw std::exception("key not found");
		}
		return value;
	}

	inline void printHierarchy(IObject o) {
		auto header = o.getHeader();
		printf("%s\n", header.getFullName().c_str());

		if (IPolyMesh::matches(header)) {
			IPolyMesh meshyObj(o);
			ICompoundProperty props = meshyObj.getProperties();
			printProperties(props);
		}

		if (ICamera::matches(header)) {
			ICamera cameraObj(o);
			ICompoundProperty props = cameraObj.getProperties();
			printProperties(props);

			ICameraSchema schema = cameraObj.getSchema();
			CameraSample sample;
			schema.get(sample);
			printf("-camera-\n");
			printf("near: %f\n", sample.getNearClippingPlane());
			printf("far: %f\n", sample.getFarClippingPlane());

			printf("focal length: %f\n", sample.getFocalLength());
			printf("horizontal aperture: %f\n", sample.getHorizontalAperture());
			printf("-camera-\n");
		}

		for (int i = 0; i < o.getNumChildren(); ++i) {
			IObject child = o.getChild(i);
			printHierarchy(child);
		}
	}

	glm::dmat4 GetTransform(IObject o) {
		Abc::M44d m;
		while (o) {
			if (IXform::matches(o.getHeader())) {
				IXform form(o);
				auto schema = form.getSchema();
				auto value = schema.getValue();
				auto matrix = value.getMatrix();
				m = m * matrix;
			}
			o = o.getParent();
		}
		glm::dmat4 matrix;
		memcpy(glm::value_ptr(matrix), m.getValue(), sizeof(double) * 16);
		return matrix;
	}

	typedef strict_variant::variant<glm::vec3, float, std::string> AttributeVariant;
	struct AlembicGeometry {
		std::vector<glm::vec3> points;
		std::vector<glm::ivec3> primitives;
		std::map<std::string, std::vector<AttributeVariant>> primitiveAttributes;
	};

	inline void parseHierarchy(IObject o, Scene &scene) {
		auto header = o.getHeader();

		if (IPolyMesh::matches(header)) {
			IPolyMesh polyMesh(o);
			IPolyMeshSchema &mesh = polyMesh.getSchema();

			auto transform = GetTransform(o);

			Geometry geometry;

			// Parse Point
			Abc::IP3fArrayProperty P = mesh.getPositionsProperty();
			P3fArraySamplePtr PSample;
			P.get(PSample);

			geometry.points.resize(PSample->size());
			for (int i = 0; i < PSample->size(); ++i) {
				auto p = PSample->get()[i];
				auto point = glm::vec3(p.x, p.y, p.z);
				point = transform * glm::vec4(point, 1.0f);
				geometry.points[i].P = glm::vec3(point.x, point.y, point.z);
			}

			Abc::IInt32ArrayProperty FaceCounts = mesh.getFaceCountsProperty();
			Int32ArraySamplePtr FaceCountsSample;
			FaceCounts.get(FaceCountsSample);

			Abc::IInt32ArrayProperty Indices = mesh.getFaceIndicesProperty();
			Int32ArraySamplePtr IndicesSample;
			Indices.get(IndicesSample);

			// Parse Attributes
			ICompoundProperty props = polyMesh.getProperties();

			std::vector<string> materials;
			std::vector<glm::vec3> Cd;
			std::vector<glm::vec3> Le;
			try {
				// Material
				materials = propertyMaterial(props);

				// Diffuse
				Cd = propertyArrayVec3(props, "/.geom/.arbGeomParams/Cd");
				Le = propertyArrayVec3(props, "/.geom/.arbGeomParams/Le");
			}
			catch (...) {
				printf("material property error\n");
				materials.resize(FaceCountsSample->size(), LambertianMaterialString);
				Cd.resize(FaceCountsSample->size(), glm::vec3(0.5f));
				Le.resize(FaceCountsSample->size(), glm::vec3());
			}


			const int32_t *indices = IndicesSample->get();

			for (int i = 0; i < FaceCountsSample->size(); ++i) {
				auto material = materials[i];
				auto count = FaceCountsSample->get()[i];

				for (int j = 2; j < count; ++j) {
					Geometry::Primitive primitive;
					primitive.indices[0] = indices[0];
					//primitive.indices[1] = indices[j - 1];
					//primitive.indices[2] = indices[j];
					primitive.indices[1] = indices[j];
					primitive.indices[2] = indices[j - 1];

					if (material == LambertianMaterialString) {
						LambertianMaterial m(Le[i], Cd[i]);
						primitive.material = m;
					}
					geometry.primitives.push_back(primitive);
				}
				indices += count;
			}

			scene.geometries.push_back(geometry);
		}

		if (ICamera::matches(header)) {
			/*
			‚Æ‚è‚ ‚¦‚¸
			@@Transform
				View/Resolution
				View/Focal Length, Aperture (fovy ‚ðŒvŽZ)
				    
				Bokeh
					Focus Distance
					F-Stop
			‚Ì‚Ý‘Î‰ž
			*/

			ICamera cameraObj(o);
			ICompoundProperty props = cameraObj.getProperties();
			ICameraSchema schema = cameraObj.getSchema();

			try {
				CameraSample sample;
				schema.get(sample);

				CameraSetting setting;
				setting.imageWidth = propertyScalarFloat(props, "/.geom/.userProperties/resx");
				setting.imageHeight = propertyScalarFloat(props, "/.geom/.userProperties/resy");

				// http://127.0.0.1:48626/nodes/obj/cam#aperture
				double aperture = sample.getVerticalAperture() /*centimeters*/ / 100.0;
				double focalLength = sample.getFocalLength() /*millimeters*/ / 1000.0;
				double fovy = std::atan2(aperture * 0.5, focalLength) * 2.0;
				setting.fovy = fovy;
				setting.focasDistance = sample.getFocusDistance();

				// FStop = FocalLength / (Radius * 2)
				// Radius = FocalLength / (2 * FStop)
				double fStop = sample.getFStop();

				printf("A, focusDistance : %f\n", setting.focasDistance);
				printf("B, focalLength : %f\n", focalLength);
				printf("f-number: %f\n", fStop);

				// Octane Setting
				// double phi = focalLength / fStop;
				// setting.lensRadius = phi;

				// Mantra setting
				{
					double A = setting.focasDistance;
					double B = focalLength;
					double f = A * B / (A + B);
					printf("f : %f\n", f);
					double phi = f / fStop;
					setting.lensRadius = phi * 0.5;
				}

				printf("setting.lensRadius : %f\n", setting.lensRadius);
				// setting.lensRadius = focusDistance / (2.0 * fStop);

				glm::dmat4 transform = GetTransform(o);
				glm::dmat4 inverseTransposed = glm::inverseTranspose(transform);

				glm::dvec4 origin = transform * glm::dvec4(0.0, 0.0, 0.0, 1.0);
				glm::dvec4 front = inverseTransposed * glm::dvec4(0.0, 0.0, -1.0, 1.0);
				glm::dvec4 up = inverseTransposed * glm::dvec4(0.0, 1.0, 0.0, 1.0);
				setting.eye = origin;
				setting.lookat = origin + front;
				setting.up = up;

				// setting.lensRadius = 0.0f;

				scene.camera = Camera(setting);
			}
			catch (...) {
				printf("camera property error\n");
			}
		}

		for (int i = 0; i < o.getNumChildren(); ++i) {
			IObject child = o.getChild(i);
			parseHierarchy(child, scene);
		}
	}

	inline void loadFromABC(const char *filename, Scene &scene) {
		IArchive archive(Alembic::AbcCoreOgawa::ReadArchive(), filename);
		IObject top(archive, kTop);
		rt::printHierarchy(top);

		rt::parseHierarchy(top, scene);
	}
}