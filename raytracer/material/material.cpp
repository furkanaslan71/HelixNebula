#include "material.h"

Material::Material()
	: id(-1),
		type(""),
		ambient_reflectance(),
		diffuse_reflectance(),
		specular_reflectance(),
		mirror_reflectance(),
		phong_exponent(0.0f),
		refraction_index(1.0f),
		absorption_coefficient(),
	absorption_index(),
	roughness()
{}

Material::Material(const Material_& _material)
		: id(_material.id),
			type(_material.type),
			ambient_reflectance(_material.ambient_reflectance),
			diffuse_reflectance(_material.diffuse_reflectance),
			specular_reflectance(_material.specular_reflectance),
			mirror_reflectance(_material.mirror_reflectance),
			phong_exponent(_material.phong_exponent),
			refraction_index(_material.refraction_index),
			absorption_coefficient(_material.absorption_coefficient),
			absorption_index(_material.absorption_index),
			roughness(_material.roughness)
{
}
