#pragma once

#include <Alembic/Abc/All.h>
#include <Alembic/AbcCoreOgawa/All.h>
#include <Alembic/AbcGeom/All.h>

#include "render_object.hpp"

#include <functional>

namespace rt {
	typedef strict_variant::variant<
		glm::dvec3,
		double,
		std::string
	> AttributeVariant;
	
	struct AlembicGeometry {
		std::vector<glm::dvec3> points;
		std::vector<glm::ivec3> primitives;
		std::map<std::string, std::vector<AttributeVariant>> primitiveAttributes;
	};

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

	inline std::vector<AttributeVariant> readAttributes(IArrayProperty prop) {
		auto dataType = prop.getDataType();
		if (dataType.getExtent() == 3 && dataType.getPod() == kFloat32POD) {
			Abc::IV3fArrayProperty arrayProp(prop.getParent(), prop.getName());
			V3fArraySamplePtr sample;
			arrayProp.get(sample);
			std::vector<AttributeVariant> values(sample->size());
			for (int j = 0; j < sample->size(); ++j) {
				auto value = sample->get()[j];
				values[j] = glm::dvec3(value.x, value.y, value.z);
			}
			return values;
		} else if(dataType.getExtent() == 1 && dataType.getPod() == kFloat32POD) {
			Abc::IFloatArrayProperty arrayProp(prop.getParent(), prop.getName());
			FloatArraySamplePtr sample;
			arrayProp.get(sample);
			std::vector<AttributeVariant> values(sample->size());
			for (int j = 0; j < sample->size(); ++j) {
				auto value = sample->get()[j];
				values[j] = value;
			}
			return values;
		}
		else {
			throw std::exception("unsupported dataType");
		}
	}
	inline std::vector<AttributeVariant> readAttributes(ICompoundProperty prop) {
		Abc::IUInt32ArrayProperty indicesProp(prop, ".indices");
		Abc::IStringArrayProperty stringProp(prop, ".vals");

		UInt32ArraySamplePtr indicesSample;
		indicesProp.get(indicesSample);
		std::vector<uint32_t> indices(indicesSample->get(), indicesSample->get() + indicesSample->size());

		StringArraySamplePtr stringSample;
		stringProp.get(stringSample);
		std::vector<std::string> values(stringSample->get(), stringSample->get() + stringSample->size());

		std::vector<AttributeVariant> expands(indices.size());
		for (int i = 0; i < indices.size(); ++i) {
			expands[i] = values[indices[i]];
		}
		return expands;
	}

	inline std::map<std::string, std::vector<AttributeVariant>> arbGeomParamsAttributes(ICompoundProperty props) {
		std::map<std::string, std::vector<AttributeVariant>> attributes;
		try {
			ICompoundProperty geom(props, ".geom");
			ICompoundProperty arbGeomParams(geom, ".arbGeomParams");
			for (int i = 0; i < arbGeomParams.getNumProperties(); ++i) {
				auto header = arbGeomParams.getPropertyHeader(i);
				auto attributeHeader = arbGeomParams.getPropertyHeader(i);
				auto propType = attributeHeader.getPropertyType();

				try {
					if (propType == AbcA::PropertyType::kScalarProperty) {
						// invalid
						continue;
					}
					else if (propType == AbcA::PropertyType::kArrayProperty) {
						auto valueProp = IArrayProperty(arbGeomParams, attributeHeader.getName());
						attributes[attributeHeader.getName()] = readAttributes(valueProp);
					}
					else if (propType == AbcA::PropertyType::kCompoundProperty) {
						auto valueProp = ICompoundProperty(arbGeomParams, attributeHeader.getName());
						attributes[attributeHeader.getName()] = readAttributes(valueProp);
					}
				}
				catch (std::exception &e) {
					printf("exception, %s\n", e.what());
					continue;
				}
			}
		}
		catch (std::exception &e) {
			printf("exception, %s\n", e.what());
		}
		return attributes;
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

	inline std::vector<glm::dvec3> toVec3Array(V3fArraySamplePtr sample) {
		std::vector<glm::dvec3> values(sample->size());
		for (int j = 0; j < sample->size(); ++j) {
			auto value = sample->get()[j];
			values[j] = glm::dvec3(value.x, value.y, value.z);
		}
		return values;
	}

	inline std::vector<glm::dvec3> propertyArrayVec3(ICompoundProperty props, const char *dir) throw(std::exception) {
		bool found = false;
		std::vector<glm::dvec3> values;
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
					values[j] = glm::dvec3(value.x, value.y, value.z);
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
	inline double propertyScalarFloat(ICompoundProperty props, const char *dir) {
		bool found = false;
		float value = 0.0f;
		visitProperties(props,
			[&](IScalarProperty prop, std::string name) {
			if (name == dir) {
				found = true;
				IFloatProperty doubleProp(prop.getParent(), prop.getName());
				doubleProp.get(value);
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

	inline void parseHierarchy(IObject o, Scene &scene, std::function<Geometry (AlembicGeometry)> binding) {
		auto header = o.getHeader();

		if (IPolyMesh::matches(header)) {
			IPolyMesh polyMesh(o);
			IPolyMeshSchema &mesh = polyMesh.getSchema();

			auto transform = GetTransform(o);
			
			AlembicGeometry geometry;

			// Parse Point
			Abc::IP3fArrayProperty P = mesh.getPositionsProperty();
			P3fArraySamplePtr PSample;
			P.get(PSample);

			geometry.points.resize(PSample->size());
			for (int i = 0; i < PSample->size(); ++i) {
				auto p = PSample->get()[i];
				auto point = glm::dvec3(p.x, p.y, p.z);
				point = transform * glm::dvec4(point, 1.0);
				geometry.points[i] = glm::dvec3(point.x, point.y, point.z);
			}

			Abc::IInt32ArrayProperty FaceCounts = mesh.getFaceCountsProperty();
			Int32ArraySamplePtr FaceCountsSample;
			FaceCounts.get(FaceCountsSample);

			Abc::IInt32ArrayProperty Indices = mesh.getFaceIndicesProperty();
			Int32ArraySamplePtr IndicesSample;
			Indices.get(IndicesSample);

			// Attributesをパース、3角形ポリゴンにする関係で、アトリビュートの配列も調整する
			ICompoundProperty props = polyMesh.getProperties();
			auto attributes = arbGeomParamsAttributes(props);
			std::map<std::string, std::vector<AttributeVariant>> primAttributes;

			const int32_t *indices = IndicesSample->get();
			//int primitiveCount = 0;
			//for (int i = 0; i < FaceCountsSample->size(); ++i) {
			//	primitiveCount += FaceCountsSample->get()[i];
			//}

			for (int i = 0; i < FaceCountsSample->size(); ++i) {
				auto count = FaceCountsSample->get()[i];

				for (int j = 2; j < count; ++j) {
					glm::ivec3 primitive;
					primitive[0] = indices[0];
					primitive[1] = indices[j];
					primitive[2] = indices[j - 1];
					geometry.primitives.push_back(primitive);

					for (auto it = attributes.begin(); it != attributes.end(); ++it) {
						if (it->second.size() == FaceCountsSample->size()) {
							AttributeVariant attrib = it->second[i];
							primAttributes[it->first].push_back(attrib);
						}
					}
				}
				indices += count;
			}
			geometry.primitiveAttributes = primAttributes;

			scene.geometries.push_back(binding(geometry));
		}

		if (ICamera::matches(header)) {
			/*
			とりあえず
			　　Transform
				View/Resolution
				View/Focal Length, Aperture (fovy を計算)
				    
				Bokeh
					Focus Distance
					F-Stop
			のみ対応
			*/

			ICamera cameraObj(o);
			ICompoundProperty props = cameraObj.getProperties();
			ICameraSchema schema = cameraObj.getSchema();

			try {
				CameraSample sample;
				schema.get(sample);

				CameraSetting setting;
				setting.imageWidth = (int)propertyScalarFloat(props, "/.geom/.userProperties/resx");
				setting.imageHeight = (int)propertyScalarFloat(props, "/.geom/.userProperties/resy");

				// http://127.0.0.1:48626/nodes/obj/cam#aperture
				double aperture = sample.getVerticalAperture() /*centimeters*/ / 100.0;
				double focalLength = sample.getFocalLength() /*millimeters*/ / 1000.0;
				double fovy = std::atan2(aperture * 0.5, focalLength) * 2.0;
				setting.fovy = (double)fovy;
				setting.focasDistance = (double)sample.getFocusDistance();

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

				// setting.lensRadius = 0.0;

				scene.camera = Camera(setting);
			}
			catch (...) {
				printf("camera property error\n");
			}
		}

		for (int i = 0; i < o.getNumChildren(); ++i) {
			IObject child = o.getChild(i);
			parseHierarchy(child, scene, binding);
		}
	}

	inline Geometry geometryMaterialBinding(AlembicGeometry abcGeom) {
		Geometry geom;
		for (int pointID = 0; pointID < abcGeom.points.size(); ++pointID) {
			Geometry::Point point;
			point.P = abcGeom.points[pointID];
			geom.points.push_back(point);
		}
		for (int primID = 0; primID < abcGeom.primitives.size(); ++primID) {
			Geometry::Primitive prim;
			prim.indices = abcGeom.primitives[primID];
			prim.material = LambertianMaterial();
			geom.primitives.push_back(prim);
		}

		// Primitive Material
		const char *LambertianMaterialString = "LambertianMaterial";
		const char *SpecularMaterialString = "SpecularMaterial";
		const char *MicrofacetConductorMaterialString = "MicrofacetConductorMaterial";
		const char *MicrofacetCoupledConductorMaterialString = "MicrofacetCoupledConductorMaterial";
		const char *MicrofacetCoupledDielectricsMaterialString = "MicrofacetCoupledDielectricsMaterial";
		const char *HeitzConductorMaterialString = "HeitzConductorMaterial";

		auto material = abcGeom.primitiveAttributes.find("Material");
		if (material == abcGeom.primitiveAttributes.end()) {
			return geom;
		}
		std::vector<AttributeVariant> materials = material->second;
		for (int primID = 0; primID < materials.size(); ++primID) {
			std::string materialString;
			if (auto m = strict_variant::get<std::string>(&materials[primID])) {
				materialString = *m;
			}
			else {
				printf("loadFromABC unknown material string[%s]\n", m->c_str());
			}

			if (materialString == LambertianMaterialString) {
				LambertianMaterial m;
				if (abcGeom.primitiveAttributes.count("Le")) {
					auto Le = abcGeom.primitiveAttributes["Le"][primID];
					if (auto LeVec3 = strict_variant::get<glm::dvec3>(&Le)) {
						m.Le = *LeVec3;
					}
				}
				if (abcGeom.primitiveAttributes.count("Cd")) {
					auto Cd = abcGeom.primitiveAttributes["Cd"][primID];
					if (auto CdVec3 = strict_variant::get<glm::dvec3>(&Cd)) {
						m.R = *CdVec3;
					}
				}
				geom.primitives[primID].material = m;
			}
			else if (materialString == MicrofacetConductorMaterialString) {
				MicrofacetConductorMaterial m;
				if (abcGeom.primitiveAttributes.count("roughness")) {
					const std::vector<AttributeVariant> &roughnesses = abcGeom.primitiveAttributes["roughness"];
					if (auto roughness = strict_variant::get<double>(&roughnesses[primID])) {
						m.alpha = (*roughness) * (*roughness);
					}
				}
				geom.primitives[primID].material = m;
			}
			else if (materialString == MicrofacetCoupledConductorMaterialString) {
				MicrofacetCoupledConductorMaterial m;
				if (abcGeom.primitiveAttributes.count("roughness")) {
					const std::vector<AttributeVariant> &roughnesses = abcGeom.primitiveAttributes["roughness"];
					if (auto roughness = strict_variant::get<double>(&roughnesses[primID])) {
						m.alpha = (*roughness) * (*roughness);
					}
				}
				geom.primitives[primID].material = m;
			}
			else if (materialString == MicrofacetCoupledDielectricsMaterialString) {
				MicrofacetCoupledDielectricsMaterial m;
				if (abcGeom.primitiveAttributes.count("roughness")) {
					const std::vector<AttributeVariant> &roughnesses = abcGeom.primitiveAttributes["roughness"];
					if (auto roughness = strict_variant::get<double>(&roughnesses[primID])) {
						m.alpha = (*roughness) * (*roughness);
					}
				}
				if (abcGeom.primitiveAttributes.count("Cd")) {
					const std::vector<AttributeVariant> &Cd = abcGeom.primitiveAttributes["Cd"];
					if (auto CdVec3 = strict_variant::get<glm::dvec3>(&Cd[primID])) {
						m.Cd = *CdVec3;
					}
				}
				geom.primitives[primID].material = m;
			}
			else if (materialString == SpecularMaterialString) {
				geom.primitives[primID].material = SpecularMaterial();
			}
			else if (materialString == HeitzConductorMaterialString) {
				double alpha = 0.5;

				if (abcGeom.primitiveAttributes.count("roughness")) {
					const std::vector<AttributeVariant> &roughnesses = abcGeom.primitiveAttributes["roughness"];
					if (auto roughness = strict_variant::get<double>(&roughnesses[primID])) {
						alpha = (*roughness) * (*roughness);
					}
				}

				geom.primitives[primID].material = HeitzConductorMaterial(alpha);
			}
		}
		return geom;
	}

	inline void loadFromABC(const char *filename, Scene &scene) {
		IArchive archive(Alembic::AbcCoreOgawa::ReadArchive(), filename);
		IObject top(archive, kTop);
		rt::printHierarchy(top);

		rt::parseHierarchy(top, scene, geometryMaterialBinding);
	}
}