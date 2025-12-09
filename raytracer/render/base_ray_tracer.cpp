#include "base_ray_tracer.h"

static double getCosTheta(glm::vec3 v1, glm::vec3 v2)
{
	v1 = glm::normalize(v1);
	v2 = glm::normalize(v2);
	return glm::dot(v1, v2);
}

static glm::vec3 reflect(glm::vec3 wo, glm::vec3 n)
{
	wo = glm::normalize(wo);

	glm::vec3 wr = (n * (2 * (glm::dot(n, wo)))) - wo;
	return wr;
}

static double getCosPhi(double cosTheta, double n1, double n2)
{
	double sinTheta2 = 1 - cosTheta * cosTheta;
	double n = n1 / n2;
	return sqrt(1 - n * n * sinTheta2);
}

static bool couldRefract(double cosTheta, double n1, double n2)
{
	double sinTheta2 = 1 - cosTheta * cosTheta;
	double n = n1 / n2;
	return (1 - n * n * sinTheta2) >= 0;
}

static glm::vec3 snellRefract(glm::vec3 wo, glm::vec3 n, double n1, double n2)
{
	double eta = n1 / n2;
	double cosTheta = std::clamp(glm::dot(wo, n), -1.0f, 1.0f);
	double sin2ThetaT = eta  * (1 - cosTheta * cosTheta);
	if (sin2ThetaT > 1.0)
		return glm::vec3(0); // total internal reflection
	double cosThetaT = sqrt(1.0 - sin2ThetaT);
	return -wo + n * (float)(eta * cosTheta - cosThetaT) * (float)eta;
}

static double r_parallel(double cosTheta, double cosPhi, double n1, double n2)
{
	return ((n2 * cosTheta) - (n1 * cosPhi)) / ((n2 * cosTheta) + (n1 * cosPhi));
}

static double r_perpendicular(double cosTheta, double cosPhi, double n1, double n2)
{
	return ((n1 * cosTheta) - (n2 * cosPhi)) / ((n1 * cosTheta) + (n2 * cosPhi));
}

static double fresnelReflectance(double r_parallel, double r_perpendicular)
{
	return (r_parallel * r_parallel + r_perpendicular * r_perpendicular) / 2.0;
}


BaseRayTracer::BaseRayTracer(Color& background_color,
	LightSources& light_sources,
	std::shared_ptr<BVH<TLASBox>> world,
	std::vector<Plane>& planes,
	MaterialManager& material_manager,
	RendererInfo& renderer_info)
	:background_color(background_color),
	 light_sources(light_sources),
		world(world),
	planes(planes),
		material_manager(material_manager),
	renderer_info(renderer_info)
{
}

Color BaseRayTracer::traceRay(const Ray& ray, const RenderContext& context) const
{
	// Placeholder implementation: return background color
	double min_t = INFINITY;
	return computeColor(ray, renderer_info.max_recursion_depth + 1, context);
}

Color BaseRayTracer::computeColor(const Ray& ray, int depth, const RenderContext& context) const
{
	if (depth <= 0) return Color(0, 0, 0);

	HitRecord rec;
	bool hit_plane = false;
	hit_plane = this->hitPlanes(ray, Interval(renderer_info.intersection_test_epsilon, INFINITY), rec);

	double closest_t = hit_plane ? rec.t : INFINITY;

	if (!world->intersect(ray, Interval(renderer_info.intersection_test_epsilon, closest_t), rec))
	{
		if (!hit_plane)
		{
			if (depth == renderer_info.max_recursion_depth + 1)
				return Color(background_color);
			else
				return Color(0, 0, 0);
		}
	}
#if BACKFACE_CULLING
	if (!ray.inside && !rec.front_face)
		return Color(0, 0, 0);
#endif
	return applyShading(ray, depth, rec, context);
}

Color BaseRayTracer::applyShading(const Ray& ray, 
	int depth, HitRecord& rec, const RenderContext& context) const
{
	Material mat = material_manager.getMaterialById(rec.material_id);
	Color color = Color(mat.ambient_reflectance) * Color(light_sources.ambient_light);

	if ((mat.type).compare("mirror") == 0)
	{
		glm::vec3 wo = -ray.direction;
		glm::vec3 wr = (rec.normal * (2 * (glm::dot(rec.normal, wo)))) - wo;
		Ray reflectedRay = Ray(rec.point + rec.normal * renderer_info.shadow_ray_epsilon, wr, ray.time);

		if (mat.roughness != 0)
		{
			reflectedRay.perturb(mat.roughness);
		}
		
		color += computeColor(reflectedRay, depth - 1, context) * Color(mat.mirror_reflectance);
	}
	else if ((mat.type).compare("conductor") == 0)
	{
		glm::vec3 wo = -ray.direction;
		glm::vec3 wr = (rec.normal * (2 * (glm::dot(rec.normal, wo)))) - wo;
		wr = glm::normalize(wr);
		wo = glm::normalize(wo);
		double cos_theta = glm::dot(wo, rec.normal);
		double k = mat.absorption_index; // Assuming k is the same for r, g, b
		double n = static_cast<double>(mat.refraction_index);
		double rs_num = (n * n) + (k * k)
			- (n * cos_theta * static_cast<double>(2.0))
			+ (cos_theta * cos_theta);
		double rs_den = (n * n) + (k * k)
			+ (n * cos_theta * static_cast<double>(2.0))
			+ (cos_theta * cos_theta);
		double rs = rs_num / rs_den;

		double rp_num = ((n * n) + (k * k)) * (cos_theta * cos_theta)
			- (n * cos_theta * static_cast<double>(2.0))
			+ 1.0;
		double rp_den = ((n * n) + (k * k)) * (cos_theta * cos_theta)
			+ (n * cos_theta * static_cast<double>(2.0))
			+ 1.0;
		double rp = rp_num / rp_den;
		double f_r = (rs + rp) * 0.5;

		Ray reflectedRay = Ray(rec.point + rec.normal * renderer_info.shadow_ray_epsilon, wr, ray.time);
		if (mat.roughness != 0)
		{
			reflectedRay.perturb(mat.roughness);
		}

		color += computeColor(reflectedRay, depth - 1, context) * f_r * mat.mirror_reflectance;
	}
	else if (mat.type == "dielectric")
	{
		glm::vec3 wo = -ray.direction;
		glm::vec3 normal = rec.normal;
		double n1, n2;
		//bool entering = rec.front_face;
		bool entering = !ray.inside;

		if (entering) { n1 = 1.0; n2 = mat.refraction_index; }
		else { n1 = mat.refraction_index; n2 = 1.0; normal = normal * -1.0f; }

		double eta = n1 / n2;
		double cosTheta = std::clamp(glm::dot(wo, normal), -1.0f, 1.0f);
		double sin2ThetaT = eta * eta * (1 - cosTheta * cosTheta);
		double F_r = 1.0;

		if (sin2ThetaT <= 1.0)
		{
			double cosThetaT = sqrt(1.0 - sin2ThetaT);
			double r_par = r_parallel(cosTheta, cosThetaT, n1, n2);
			double r_perp = r_perpendicular(cosTheta, cosThetaT, n1, n2);
			F_r = fresnelReflectance(r_par, r_perp);
		}

		// Reflection
		glm::vec3 wr = glm::normalize(reflect(wo, normal));

		Ray reflectedRay(rec.point + normal * renderer_info.shadow_ray_epsilon, wr, ray.time);
		reflectedRay.inside = ray.inside;

		if (mat.roughness != 0)
		{
			reflectedRay.perturb(mat.roughness);
		}

		// Refraction
		Color refractedColor(0.0, 0.0, 0.0);
		if (sin2ThetaT <= 1.0)
		{
			glm::vec3 wt = -wo * (float)eta + normal * (float)(eta * cosTheta - sqrt(1 - sin2ThetaT));
			wt = glm::normalize(wt);
			Ray refractedRay(rec.point - normal * renderer_info.shadow_ray_epsilon, wt, ray.time);
			refractedRay.inside = !ray.inside;
			if (mat.roughness != 0)
			{
				refractedRay.perturb(mat.roughness);
			}
			refractedColor = computeColor(refractedRay, depth - 1, context);
		}

		Color reflectedColor = computeColor(reflectedRay, depth - 1, context);
		//reflectedColor = Color(0.0);
		Color L = reflectedColor * F_r
			+ refractedColor * (1 - F_r)
			;


		// Absorption when exiting
		if (!entering)
		{
			double d = rec.t; // or track actual thickness
			L.r *= exp(-mat.absorption_coefficient.x * d);
			L.g *= exp(-mat.absorption_coefficient.y * d);
			L.b *= exp(-mat.absorption_coefficient.z * d);
		}
		//L = Color(0.0);
		return color + L;
	}

	for (const auto& light : light_sources.point_lights)
	{
		glm::vec3 wi = glm::vec3(light.position) - rec.point;
		double distance = glm::length(wi);
		wi = glm::normalize(wi);
		Ray shadowRay = Ray(rec.point + rec.normal * renderer_info.shadow_ray_epsilon, wi, ray.time);
		Ray shadowRayPlane = Ray(rec.point 
			+ rec.normal 
			* static_cast<float>(renderer_info.shadow_ray_epsilon) 
			* 1e-2f
			, wi, ray.time);
		HitRecord shadowRec;
		HitRecord planeShadowRec;
		if (!world->intersect(shadowRay, Interval(renderer_info.intersection_test_epsilon, distance), shadowRec) 
			&& !hitPlanes(shadowRayPlane, Interval(0, distance), planeShadowRec))
		{
			// Diffuse
			double cosTheta = std::max(0.0f, glm::dot(rec.normal, wi));
			color += Color(mat.diffuse_reflectance) * Color(light.intensity) * (cosTheta / (distance * distance));

			// Specular
			glm::vec3 wo = (ray.origin - rec.point);
			wo = glm::normalize(wo);
			glm::vec3 h = (wi + wo);
			h = glm::normalize(h);
			double cosAlpha = std::max(0.0f, glm::dot(rec.normal, h));
			color += Color(mat.specular_reflectance) * Color(light.intensity) * (pow(cosAlpha, mat.phong_exponent) / (distance * distance));
		}
	}
	for (const auto& light : light_sources.area_lights)
	{
		glm::vec3 light_sample_point = (*context.area_light_samples)[depth - 1][light.id][context.sample_index];

		// Calculate Shadow Ray towards this specific random point
		glm::vec3 wi = light_sample_point - rec.point;
		double distance = glm::length(wi);
		wi = glm::normalize(wi);
		Ray shadowRay = Ray(rec.point + rec.normal * renderer_info.shadow_ray_epsilon, wi, ray.time);
		Ray shadowRayPlane = Ray(rec.point
			+ rec.normal
			* static_cast<float>(renderer_info.shadow_ray_epsilon)
			* 1e-2f
			, wi, ray.time);
		HitRecord shadowRec;
		HitRecord planeShadowRec;
		if (!world->intersect(shadowRay, Interval(renderer_info.intersection_test_epsilon, distance), shadowRec)
			&& !hitPlanes(shadowRayPlane, Interval(0, distance), planeShadowRec))
		{
			float area_of_light = light.edge * light.edge;
			double cos_alpha_light = std::abs(glm::dot(-wi, light.normal));
			float distance2 = distance * distance;
			float attenuation = area_of_light * cos_alpha_light / distance2;
			// Diffuse
			double cosTheta = std::max(0.0f, glm::dot(rec.normal, wi));
			color += Color(mat.diffuse_reflectance) * Color(light.radiance) * attenuation * cosTheta;

			// Specular
			glm::vec3 wo = (ray.origin - rec.point);
			wo = glm::normalize(wo);
			glm::vec3 h = (wi + wo);
			h = glm::normalize(h);
			double cosAlpha = std::max(0.0f, glm::dot(rec.normal, h));
			color += Color(mat.specular_reflectance) * Color(light.radiance) *
				attenuation * pow(cosAlpha, mat.phong_exponent);
		}
	}
	return color;
}

bool BaseRayTracer::hitPlanes(const Ray& ray, Interval ray_t, HitRecord& rec) const
{
	bool hit_anything = false;
	auto closest_so_far = ray_t.max;
	HitRecord temp_rec;

	for (const auto& plane : planes)
	{
		if (plane.hit(ray, Interval(ray_t.min, closest_so_far), temp_rec))
		{
			hit_anything = true;
			closest_so_far = temp_rec.t; 
			rec = temp_rec;              
		}
	}

	return hit_anything;
}



