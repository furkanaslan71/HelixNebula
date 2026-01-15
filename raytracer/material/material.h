#ifndef MATERIAL_H
#define MATERIAL_H
#include <string>

#include "parser/parser.hpp"

enum class BrdfType : int {
	OriginalBlinnPhong,
	OriginalPhong,
	ModifiedBlinnPhong,
	ModifiedPhong,
	TorranceSparrow
};

struct BRDF {
	BRDF() = default;
	BRDF(const BRDF_& _brdf);
	BrdfType type;
	int id;
	bool normalized;
	float exponent;
	bool kdfresnel;
	glm::vec3 Evaluate(const glm::vec3& wi, const glm::vec3& wo, const glm::vec3& n,
					   const glm::vec3& kd, const glm::vec3& ks) const;
};

class Material {
public:
	Material();
	Material(const Material_& _material, BRDF* brdf = nullptr);
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
	bool degamma;
	BRDF* brdf;
};


#endif // MATERIAL_H
