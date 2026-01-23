#include "material.h"

BRDF::BRDF(const BRDF_& _brdf)
	: id(_brdf.id), exponent(_brdf.exponent), kdfresnel(_brdf.kdfresnel), normalized(_brdf.normalized)
{
	if (_brdf.type == "OriginalBlinnPhong")
	{
		type = BrdfType::OriginalBlinnPhong;
	}
	else if (_brdf.type == "OriginalPhong")
	{
		type = BrdfType::OriginalPhong;
	}
	else if (_brdf.type == "ModifiedBlinnPhong")
	{
		type = BrdfType::ModifiedBlinnPhong;
	}
	else if (_brdf.type == "ModifiedPhong")
	{
		type = BrdfType::ModifiedPhong;
	}
	else if (_brdf.type == "TorranceSparrow")
	{
		type = BrdfType::TorranceSparrow;
	}
	else
	{
		throw std::runtime_error("Unsupported BRDF type!");
	}
}

glm::vec3 BRDF::Evaluate(const glm::vec3& wi, const glm::vec3& wo, const glm::vec3& n,
                         const glm::vec3& kd, const glm::vec3& ks, float refraction) const
{
    float cosThetaI = std::max(0.0f, glm::dot(n, wi));

    if (cosThetaI <= 0.0f) return glm::vec3(0.0f);

    glm::vec3 diffuse(0.0f);
    glm::vec3 specular(0.0f);

    glm::vec3 wh = glm::normalize(wi + wo);
    float cosAlphaH = std::max(0.0f, glm::dot(n, wh));

    switch (type)
    {
	    case BrdfType::OriginalPhong:
	    {
	        glm::vec3 r = glm::normalize(2.0f * glm::dot(n, wi) * n - wi);
	        float cosAlphaR = std::max(0.0f, glm::dot(r, wo));

	        diffuse = kd;
	        float spec_factor = std::pow(cosAlphaR, exponent);
	        if (cosThetaI > 1e-6f) spec_factor /= cosThetaI;

	        specular = ks * spec_factor;
	        break;
	    }
	    case BrdfType::ModifiedPhong:
	    {
	        glm::vec3 r = glm::normalize(2.0f * glm::dot(n, wi) * n - wi);
	        float cosAlphaR = std::max(0.0f, glm::dot(r, wo));

	        if (normalized)
	        {
	            diffuse = kd * (float)(1.0 / M_PI);
	            float norm_factor = (exponent + 2.0f) / (2.0f * M_PI);
	            specular = ks * (norm_factor * std::pow(cosAlphaR, exponent));
	        }
	    	else
	    	{
	            diffuse = kd;
	            specular = ks * std::pow(cosAlphaR, exponent);
	        }
	        break;
	    }
	    case BrdfType::OriginalBlinnPhong:
	    {
	        diffuse = kd;
	        float spec_factor = std::pow(cosAlphaH, exponent);
	        if (cosThetaI > 1e-6f) spec_factor /= cosThetaI;

	        specular = ks * spec_factor;
	        break;
	    }
	    case BrdfType::ModifiedBlinnPhong:
	    {
	        if (normalized)
	        {
	            diffuse = kd * (float)(1.0 / M_PI);
	            float norm_factor = (exponent + 8.0f) / (8.0f * M_PI);
	            specular = ks * (norm_factor * std::pow(cosAlphaH, exponent));
	        }
	    	else
	    	{
	            diffuse = kd;
	            specular = ks * std::pow(cosAlphaH, exponent);
	        }
	        break;
	    }
	    case BrdfType::TorranceSparrow:
	    {
    		float D = ((exponent + 2.0f) / (2.0f * M_PI)) * std::pow(cosAlphaH, exponent);
    		float cosBeta = std::max(0.0f, glm::dot(wi, wh));
    		float r0_val = std::pow(refraction - 1.0f, 2.0f) / std::pow(refraction + 1.0f, 2.0f);

    		glm::vec3 R0 = ks * r0_val;
    		glm::vec3 F = R0 + (glm::vec3(1.0f) - R0) * std::pow(1.0f - cosBeta, 5.0f);

    		float NdotWh = cosAlphaH;
    		float NdotWo = std::max(0.0f, glm::dot(n, wo));
    		float NdotWi = cosThetaI; // Already calculated at top
    		float WoDotWh = std::max(1e-6f, glm::dot(wo, wh));

    		float G = std::min(1.0f, std::min((2.0f * NdotWh * NdotWo) / WoDotWh,
											 (2.0f * NdotWh * NdotWi) / WoDotWh));

    		if (NdotWo > 1e-6f && NdotWi > 1e-6f)
    		{
    			specular = (D * G * F) / (4.0f * NdotWo * NdotWi);
    		}
    		else
    		{
    			specular = glm::vec3(0.0f);
    		}

    		diffuse = kd * (float)(1.0 / M_PI);

    		if (kdfresnel)
    		{
    			diffuse *= (glm::vec3(1.0f) - F);
    		}
    		break;
	    }
    }

    return diffuse + specular;
}

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
	roughness(),
	degamma(),
	scattering_coefficient(),
	anisotropy()
{}

Material::Material(const Material_& _material, BRDF* brdf)
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
			roughness(_material.roughness),
			degamma(_material.degamma),
			brdf(brdf),
			scattering_coefficient(_material.scattering_coefficient),
			anisotropy(_material.anisotropy)
{
}
