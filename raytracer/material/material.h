#ifndef MATERIAL_H
#define MATERIAL_H

#include "../include/hittable.h"
#include "../include/ray.h"
#include "../include/color.h"
#include "../include/parser.hpp"
#include <string>
#include "../include/vec3.h"

class Material {
public:
	Material();
	Material(const Material_& _material);
	int id;
	std::string type;
	Vec3 ambient_reflectance;
	Vec3 diffuse_reflectance;
	Vec3 specular_reflectance;
	Vec3 mirror_reflectance;
	float phong_exponent;
	float refraction_index;
	Vec3 absorption_coefficient;
	float absorption_index;
	float roughness;
};


#endif // MATERIAL_H
