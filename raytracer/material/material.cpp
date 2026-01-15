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
                         const glm::vec3& kd, const glm::vec3& ks) const
{
    // Cosine of the incoming angle (light vs normal)
    float cosThetaI = std::max(0.0f, glm::dot(n, wi));

    // If light is below horizon, return black (strict check from PDF Eq 2,3,5...)
    if (cosThetaI <= 0.0f) return glm::vec3(0.0f);

    glm::vec3 diffuse(0.0f);
    glm::vec3 specular(0.0f);

    // Common term: Half Vector
    glm::vec3 wh = glm::normalize(wi + wo);
    float cosAlphaH = std::max(0.0f, glm::dot(n, wh));

    switch (type)
    {
    case BrdfType::OriginalPhong:
    {
        // Eq 2: kd + ks * (cos^p(alpha_r) / cosThetaI)
        // Reflection of wi around n: R = 2(n.wi)n - wi
        glm::vec3 r = glm::normalize(2.0f * glm::dot(n, wi) * n - wi);
        float cosAlphaR = std::max(0.0f, glm::dot(r, wo));

        diffuse = kd;
        float spec_factor = std::pow(cosAlphaR, exponent);
        // Avoid division by zero
        if (cosThetaI > 1e-6f) spec_factor /= cosThetaI;

        specular = ks * spec_factor;
        break;
    }
    case BrdfType::ModifiedPhong:
    {
        // Eq 3: kd + ks * cos^p(alpha_r) (No divide by cosThetaI)
        glm::vec3 r = glm::normalize(2.0f * glm::dot(n, wi) * n - wi);
        float cosAlphaR = std::max(0.0f, glm::dot(r, wo));

        if (normalized) {
            // Eq 5: kd/PI + ks * (p+2)/2PI * cos^n(alpha_r)
            diffuse = kd * (float)(1.0 / M_PI);
            float norm_factor = (exponent + 2.0f) / (2.0f * M_PI);
            specular = ks * (norm_factor * std::pow(cosAlphaR, exponent));
        } else {
            // Eq 3
            diffuse = kd;
            specular = ks * std::pow(cosAlphaR, exponent);
        }
        break;
    }
    case BrdfType::OriginalBlinnPhong:
    {
        // Eq 7: kd + ks * (cos^p(alpha_h) / cosThetaI)
        diffuse = kd;
        float spec_factor = std::pow(cosAlphaH, exponent);
        if (cosThetaI > 1e-6f) spec_factor /= cosThetaI;

        specular = ks * spec_factor;
        break;
    }
    case BrdfType::ModifiedBlinnPhong:
    {
        if (normalized) {
            // Eq 9: kd/PI + ks * (p+8)/8PI * cos^n(alpha_h)
            diffuse = kd * (float)(1.0 / M_PI);
            float norm_factor = (exponent + 8.0f) / (8.0f * M_PI);
            specular = ks * (norm_factor * std::pow(cosAlphaH, exponent));
        } else {
            // Eq 8: kd + ks * cos^p(alpha_h)
            diffuse = kd;
            specular = ks * std::pow(cosAlphaH, exponent);
        }
        break;
    }
    case BrdfType::TorranceSparrow:
    {
        // Eq 10
        // 1. D(alpha) - Blinn's Distribution (Eq 11)
        float D = ((exponent + 2.0f) / (2.0f * M_PI)) * std::pow(cosAlphaH, exponent);

        // 2. F(beta) - Fresnel Schlick (Eq 13, 14)
        // beta is angle between wi and wh (or wo and wh, they are symmetric)
        // However, in standard Schlick, it's angle between view and half vector.
        float cosBeta = std::max(0.0f, glm::dot(wo, wh));
        // We assume R0 comes from material IOR, but PDF says use simplified:
        // We can't access IOR here easily unless passed, assumingks is F0 or similar?
        // Usually for TS, ks represents the Fresnel response at normal incidence if explicit.
        // Let's assume ks holds the base specular color/reflectance.
        // Or calculate R0 from a default IOR (e.g. 1.5) if not provided.
        // For this HW, often R0 is derived from IOR or ks is treated as reflectance.
        // Let's implement R0 based on a generic IOR of 1.5 if not available,
        // OR better: Assume ks *is* the specular reflectance color which acts as R0 approx.
        // Implementation of Eq 13 using ks as R0 equivalent color:
        glm::vec3 F = ks + (glm::vec3(1.0f) - ks) * std::pow(1.0f - cosBeta, 5.0f);

        // 3. G(wi, wo) - Masking/Shadowing (Eq 12)
        float NdotWh = std::max(0.0f, glm::dot(n, wh));
        float NdotWo = std::max(0.0f, glm::dot(n, wo));
        float NdotWi = std::max(0.0f, glm::dot(n, wi));
        float WoDotWh = std::max(1e-6f, glm::dot(wo, wh));

        float term1 = (2.0f * NdotWh * NdotWo) / WoDotWh;
        float term2 = (2.0f * NdotWh * NdotWi) / WoDotWh;
        float G = std::min(1.0f, std::min(term1, term2));

        // Denominator: 4 * cosThetaO * cosThetaI
        float denominator = 4.0f * NdotWo * NdotWi;

        if (denominator < 1e-6f) specular = glm::vec3(0.0f);
        else specular = (D * G * F) / denominator;

        // Diffuse Term handling kdfresnel
        if (kdfresnel) {
            // (1 - F) * kd / PI
            diffuse = (glm::vec3(1.0f) - F) * kd * (float)(1.0 / M_PI);
        } else {
            diffuse = kd * (float)(1.0 / M_PI);
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
	degamma()
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
			brdf(brdf)
{
}
