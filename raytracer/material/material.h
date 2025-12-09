#ifndef MATERIAL_H
#define MATERIAL_H
#include <string>

#include "core/hittable.h"
#include "core/ray.h"
#include "core/color.h"
#include "parser/parser.hpp"
#include "math_core/math_core.h"

class Material {
public:
	Material();
	Material(const Material_& _material);
	int id;
	std::string type;
	glm::vec3 ambient_reflectance;
	glm::vec3 diffuse_reflectance;
	glm::vec3 specular_reflectance;
	glm::vec3 mirror_reflectance;
	float phong_exponent;
	float refraction_index;
	glm::vec3 absorption_coefficient;
	float absorption_index;
	float roughness;
};


#endif // MATERIAL_H
