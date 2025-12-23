#include "base_ray_tracer.h"


inline glm::vec3 get_cube_tbn_normal(
	const glm::vec3& hit_normal,
	const glm::vec3& rgb_sample,
	const glm::vec3& tangent_u,
	const glm::vec3& tangent_v)
{
	glm::vec3 n_ts = rgb_sample / 127.5f - glm::vec3(1.0f);
	n_ts = glm::normalize(n_ts);

	glm::vec3 T = glm::normalize(tangent_u);
	glm::vec3 B = glm::normalize(tangent_v);
	glm::vec3 res = { T * n_ts.x +
	B * n_ts.y +
	hit_normal * n_ts.z };
	return glm::normalize(res);
	
}

inline glm::vec3 get_sphere_tbn_normal(const glm::vec3& hit_normal,
																			 const glm::vec3& rgb_sample,
																			 const glm::vec3& tangent_u,
																			 const glm::vec3& tangent_v)
{
	glm::vec3 tangent_normal = (rgb_sample / 127.5f) - glm::vec3(1.0f);
	tangent_normal = glm::normalize(tangent_normal);
	glm::vec3 T, B;

	T = tangent_u;
	B = tangent_v;

	glm::vec3 res = (T * tangent_normal.x) +
		(B * tangent_normal.y) +
		(hit_normal * tangent_normal.z);
	res = glm::normalize(res);
	return res;
}

static inline glm::vec3 reflect(glm::vec3 wo, glm::vec3 n)
{
	//wo = glm::normalize(wo);

	glm::vec3 wr = (n * (2 * (glm::dot(n, wo)))) - wo;
	return wr;
}

static inline double r_parallel(double cosTheta, double cosPhi, double n1, double n2)
{
	return ((n2 * cosTheta) - (n1 * cosPhi)) / ((n2 * cosTheta) + (n1 * cosPhi));
}

static inline double r_perpendicular(double cosTheta, double cosPhi, double n1, double n2)
{
	return ((n1 * cosTheta) - (n2 * cosPhi)) / ((n1 * cosTheta) + (n2 * cosPhi));
}

static inline double fresnelReflectance(double r_parallel, double r_perpendicular)
{
	return (r_parallel * r_parallel + r_perpendicular * r_perpendicular) / 2.0;
}


BaseRayTracer::BaseRayTracer(Color& background_color,
	LightSources& light_sources,
	std::shared_ptr<BVH<TLASBox>> world,
	std::vector<Plane>& planes,
	MaterialManager& material_manager,
	RendererInfo& renderer_info,
	TextureFetcher& texture_fetcher)
	:background_color(background_color),
	 light_sources(light_sources),
		world(world),
	planes(planes),
		material_manager(material_manager),
	renderer_info(renderer_info),
	texture_fetcher(texture_fetcher)
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
				if (renderer_info.background_tex_id != -1)
				{
					glm::vec3 d = glm::normalize(ray.direction);
					auto lookup = texture_fetcher.get_background_texture(d, renderer_info.background_tex_id, this->bgc_data);
					return Color(lookup);
				}
				else
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
	Color color(0.0, 0.0, 0.0);

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
		//wo = glm::normalize(wo);
		double cos_theta = glm::dot(wo, rec.normal);
		double cos2 = cos_theta * cos_theta;

		double k = mat.absorption_index; // Assuming k is the same for r, g, b
		double k2 = k * k;

		double n = static_cast<double>(mat.refraction_index);
		double n2 = n * n;

		double two_n_cos = 2.0 * n * cos_theta;

		double n2_k2 = n2 + k2;

		double rs_num = n2_k2 - two_n_cos + cos2;
		double rs_den = n2_k2 + two_n_cos + cos2;
		double rs = rs_num / rs_den;

		double rp_num = n2_k2 * cos2 - two_n_cos + 1.0;
		double rp_den = n2_k2 * cos2 + two_n_cos + 1.0;
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
		glm::vec3 geometric_normal = rec.normal;

		bool entering = glm::dot(ray.direction, geometric_normal) < 0;

		double n1, n2;
		glm::vec3 normal;
		double current_ior = (double)mat.refraction_index;

		if (entering)
		{
			n1 = 1.0;
			n2 = current_ior;
			normal = geometric_normal;
		}
		else
		{
			n1 = current_ior;
			n2 = 1.0;
			normal = -geometric_normal;
		}

		double eta = n1 / n2;
		double cosTheta = std::clamp(glm::dot(wo, normal), -1.0f, 1.0f);
		double sin2ThetaT = eta * eta * (1.0 - cosTheta * cosTheta);

		double F_r = 1.0;

		bool can_refract = sin2ThetaT <= 1.0;

		if (can_refract)
		{
			double cosThetaT = std::sqrt(std::max(0.0, 1.0 - sin2ThetaT));

			double r_par = r_parallel(cosTheta, cosThetaT, n1, n2);
			double r_perp = r_perpendicular(cosTheta, cosThetaT, n1, n2);
			F_r = fresnelReflectance(r_par, r_perp);
		}

		glm::vec3 wr = reflect(wo, normal);
		Ray reflectedRay(rec.point + normal * (float)renderer_info.shadow_ray_epsilon, wr, ray.time);

		reflectedRay.inside = ray.inside;

		if (mat.roughness != 0) reflectedRay.perturb(mat.roughness);

		Color refractedColor(0, 0, 0);
		if (can_refract)
		{
			double cosThetaT = std::sqrt(std::max(0.0, 1.0 - sin2ThetaT));
			glm::vec3 wt = -wo * (float)eta + normal * (float)(eta * cosTheta - cosThetaT);
			wt = glm::normalize(wt);

			Ray refractedRay(rec.point - normal * (float)renderer_info.shadow_ray_epsilon, wt, ray.time);

			refractedRay.inside = entering;

			if (mat.roughness != 0) refractedRay.perturb(mat.roughness);

			refractedColor = computeColor(refractedRay, depth - 1, context);
		}

		Color reflectedColor = computeColor(reflectedRay, depth - 1, context);
		Color L = reflectedColor * F_r + refractedColor * (1.0 - F_r);

		if (ray.inside || !entering)
		{
			double d = rec.t;
			L.r *= std::exp(-mat.absorption_coefficient.x * d);
			L.g *= std::exp(-mat.absorption_coefficient.y * d);
			L.b *= std::exp(-mat.absorption_coefficient.z * d);
		}

		return color + L;
	}

	glm::vec3 kd = mat.diffuse_reflectance;
	glm::vec3 ks = mat.specular_reflectance;
	glm::vec3 ka = mat.ambient_reflectance;
	if (!rec.texture_ids.empty())
	{
		for (const auto& id : rec.texture_ids)
		{
			auto info = texture_fetcher.get_lookup_info(rec.uv, id, rec.point);
			Expects(info.tex_val.x + 1e-8f > 0.f || info.tex_val.x < 1.f + 1e-8f);
			Expects(info.tex_val.y + 1e-8f > 0.f || info.tex_val.y < 1.f + 1e-8f);
			Expects(info.tex_val.z + 1e-8f > 0.f || info.tex_val.z < 1.f + 1e-8f);
			switch (info.mode)
			{
			case DecalMode::replace_kd:
			{
				kd = info.tex_val;
				break;
			}
			case DecalMode::replace_ks:
			{
				ks = info.tex_val;
				break;
			}
			case DecalMode::blend_kd:
			{
				kd = (info.tex_val + kd) * 0.5f;
				break;
			}
			case DecalMode::replace_all:
			{
				return Color(info.tex_val * 255.0f);
			}
			case DecalMode::replace_normal:
			{
				glm::vec3 rgb = info.tex_val * 255.0f;
				if (rec.sphere_r != -1)
					rec.normal = get_sphere_tbn_normal(rec.normal, rgb, rec.tangent_u, rec.tangent_v);
				else
					rec.normal = get_cube_tbn_normal(rec.normal, rgb, rec.tangent_u, rec.tangent_v);
				break;
			}
			case DecalMode::bump_normal:
			{
				if (texture_fetcher.getTexType(id) == TextureType::image)
				{
					glm::vec3 n_uv = rec.normal;
					glm::vec3 p_u = rec.tangent_u;
					glm::vec3 p_v = rec.tangent_v;
					float h_u = texture_fetcher.height_function(rec.uv, id, FetchMode::derivative_u);
					float h_v = texture_fetcher.height_function(rec.uv, id, FetchMode::derivative_v);
					glm::vec3 grad = h_u * p_u + h_v * p_v;
					rec.normal = glm::normalize(n_uv - grad);
				}
				else
				{
					//perlin bump
					glm::vec3 g = texture_fetcher.getPerlinGradient(rec.point, id);

					// Gradyanýn normale dik olan bileþenini bul (Surface Gradient)
					glm::vec3 g_parallel = glm::dot(rec.normal, g) * rec.normal;
					glm::vec3 g_perp = g - g_parallel;

					// Önemli: BumpFactor ile çarparak þiddeti ayarla
					rec.normal = glm::normalize(rec.normal - (g_perp * texture_fetcher.data_.textures[id].bump_factor));
				}
				
				break;
			}
			default:
				break;

			}

		}
	}

	color += Color(ka)* Color(light_sources.ambient_light);

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
			color += Color(kd) * Color(light.intensity) * (cosTheta / (distance * distance));

			// Specular
			glm::vec3 wo = (ray.origin - rec.point);
			wo = glm::normalize(wo);
			glm::vec3 h = (wi + wo);
			h = glm::normalize(h);
			double cosAlpha = std::max(0.0f, glm::dot(rec.normal, h));
			color += Color(ks) * Color(light.intensity) * (pow(cosAlpha, mat.phong_exponent) / (distance * distance));
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
			color += Color(kd) * Color(light.radiance) * attenuation * cosTheta;

			// Specular
			glm::vec3 wo = (ray.origin - rec.point);
			wo = glm::normalize(wo);
			glm::vec3 h = (wi + wo);
			h = glm::normalize(h);
			double cosAlpha = std::max(0.0f, glm::dot(rec.normal, h));
			color += Color(ks) * Color(light.radiance) *
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



