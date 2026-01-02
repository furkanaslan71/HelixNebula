#include "raytracer.h"
#include <typeinfo>


Raytracer::Raytracer(std::unique_ptr<Scene> _scene, const RenderContext& _render_context)
    : render_context(_render_context)
{
    scene = std::move(_scene);
}

void Raytracer::renderScene() const
{
    std::string saveDir = FS::absolute(__FILE__).parent_path() / "../../outputs";
    for (const auto& cam : scene->cameras)
    {
        std::vector<std::vector<Color>> image;
        renderOneCamera(cam, image);
    	std::string image_name = cam->image_name;
    	std::string extension = getFileExtension(image_name);

    	if (extension== "png" || extension == "jpg" || extension == ".jpeg")
			saveImage(saveDir, image_name, image, ImageType::SDR);
    	else if (extension == "exr" || extension == "hdr")
    	{
    		for (const Tonemap& tm : cam->tonemaps)
    		{
    			//tonemap
    			std::vector<std::vector<Color>> tonemapped_image;
    			tonemapped_image.resize(cam->image_height, std::vector<Color>(cam->image_width, Color(0, 0, 0)));
    			tonemap(image, tonemapped_image, tm);
    			saveImage(saveDir, image_name + tm.extension, tonemapped_image, ImageType::HDR);
    		}
    	}
    	else
    		throw std::runtime_error("Unsupported image type");
    }
}

void Raytracer::renderOneCamera(std::shared_ptr<BaseCamera> camera, std::vector<std::vector<Color>>& output) const
{
    output.resize(camera->image_height, std::vector<Color>(camera->image_width, Color(0, 0, 0)));

#if !MULTI_THREADING
    // SINGLE THREADED VERSION
    renderLoop(0, 1, camera, output);
#else
    // MULTI THREADED VERSION
    const int numThreads = std::thread::hardware_concurrency();
    std::vector<std::thread> threads(numThreads);
    for (int threadId = 0; threadId < numThreads; threadId++)
    {
        threads[threadId] = std::thread(&Raytracer::renderLoop, this, threadId, numThreads, camera, std::ref(output));
    }
    for (auto& t : threads) t.join();
#endif
}

void generateRaySamples(std::shared_ptr<BaseCamera> camera, int i, int j,
                        std::vector<std::pair<glm::vec3, glm::vec3>>& out_samples,
                        std::mt19937& rng)
{
    out_samples.clear();
    auto pixel_samples = generateJitteredSamples(camera->num_samples);

    // 1. Try to cast to DistributionCamera
    if (auto distCam = std::dynamic_pointer_cast<DistributionCamera>(camera))
    {
        // Now 'distCam' is recognized as a DistributionCamera*
        auto aperture_samples = generateJitteredSamples(distCam->num_samples);
        std::shuffle(aperture_samples.begin(), aperture_samples.end(), rng);

        for (size_t k = 0; k < distCam->num_samples; k++)
        {
            // Now you can access .q, .su, .sv, .aperture_size, etc.
            glm::vec3 pixel_sample = distCam->q + distCam->su * (j + pixel_samples[k].first)
                                     + distCam->sv * (i + pixel_samples[k].second);

            glm::vec3 aperture_sample = distCam->position + (distCam->u * (aperture_samples[k].first - 0.5f)
                                        + distCam->v * (aperture_samples[k].second - 0.5f)) * distCam->aperture_size;

            glm::vec3 a = (distCam->aperture_size > 0.0) ? aperture_sample : distCam->position;

            glm::vec3 dir = pixel_sample - distCam->position;
            float t_fp = distCam->focus_distance / glm::dot(dir, -distCam->w);
            glm::vec3 p = distCam->position + dir * t_fp;
            glm::vec3 d = glm::normalize(p - a);

            out_samples.emplace_back(a, d);
        }
    }
    // 2. Try to cast to PinholeCamera
    else if (auto pinCam = std::dynamic_pointer_cast<PinholeCamera>(camera))
    {
        for (size_t k = 0; k < pinCam->num_samples; k++)
        {
            glm::vec3 pixel_sample = pinCam->q + pinCam->su * (j + pixel_samples[k].first)
                                     + pinCam->sv * (i + pixel_samples[k].second);
            glm::vec3 dir = glm::normalize(pixel_sample - pinCam->position);
            out_samples.emplace_back(pinCam->position, dir); // Usually pinhole starts at position
        }
    }
}


// Helper to keep the actual rendering logic in one place
void Raytracer::renderLoop(int threadId, int stride, std::shared_ptr<BaseCamera> camera, std::vector<std::vector<Color>>& output) const
{
    std::vector<std::pair<glm::vec3, glm::vec3>> ray_samples;
    ray_samples.reserve(camera->num_samples);

    // Initialize area light sample buffers
    size_t area_light_size = scene->light_sources.area_lights.size();
    std::vector<std::vector<std::vector<glm::vec3>>> area_light_samples(render_context.max_recursion_depth + 1);
    for (auto& depth_layer : area_light_samples) {
        depth_layer.resize(area_light_size);
        for (auto& light_buffer : depth_layer) light_buffer.reserve(camera->num_samples);
    }

    std::mt19937 rng(std::random_device{}());

    for (int i = threadId; i < camera->image_height; i += stride)
    {
        for (int j = 0; j < camera->image_width; ++j)
        {
#if DEBUG_PIXEL
            if (i != HEIGHT || j != WIDTH) continue; // Only process the debug pixel
#endif
            // 1. Generate Light Samples for this pixel
            for (auto& depth : area_light_samples) {
                for (int f = 0; f < area_light_size; f++) {
                    scene->light_sources.area_lights[f].generateAreaLightSamples(depth[f], camera->num_samples);
                    std::shuffle(depth[f].begin(), depth[f].end(), rng);
                }
            }

            // 2. Generate Primary Ray Samples (Polymorphic)
            generateRaySamples(camera, i, j, ray_samples, rng);

            // 3. Trace Rays
            Color pixel_color(0.0, 0.0, 0.0);
            SamplingContext sampling_context;
            sampling_context.area_light_samples = &area_light_samples;

            for (int k = 0; k < camera->num_samples; k++)
            {
                sampling_context.sample_index = k;
                Ray primary_ray(ray_samples[k].first, ray_samples[k].second, generateRandomFloat(0, 1));
                pixel_color += traceRay(primary_ray, sampling_context, camera->context);
            }

            output[i][j] = (pixel_color / (float)camera->num_samples);
        }
    }
}
static inline glm::vec3 reflect(glm::vec3 wo, glm::vec3 n)
{
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

Color Raytracer::traceRay(const Ray& ray, const SamplingContext& sampling_context, const CameraContext& cam_context) const
{
	// Placeholder implementation: return background color
	double min_t = INFINITY;
	return computeColor(ray, render_context.max_recursion_depth + 1, sampling_context, cam_context);
}

Color Raytracer::computeColor(const Ray& ray, int depth, const SamplingContext& sampling_context, const CameraContext& cam_context) const
{
	if (depth <= 0) return Color(0, 0, 0);

	HitRecord rec;
	bool hit_plane = false;
	hit_plane = this->hitPlanes(ray, Interval(render_context.intersection_test_epsilon, INFINITY), rec);

	double closest_t = hit_plane ? rec.t : INFINITY;

	if (!scene->world->intersect<false>(ray, Interval(render_context.intersection_test_epsilon, closest_t), rec))
	{
		if (!hit_plane)
		{
			if (depth == render_context.max_recursion_depth + 1)
			{
				if (render_context.b_type == BackgroundType::Color)
					return Color(render_context.background_info.background_color);

				return Color(lookupBackgroundTex(render_context.background_info.background_tex,
					glm::normalize(ray.direction), cam_context));

			}
			return Color(0, 0, 0);
		}
	}
#if BACKFACE_CULLING
	if (!ray.inside && !rec.front_face)
		return Color(0, 0, 0);
#endif

	return applyShading(ray, depth, rec, sampling_context, cam_context);
}

Color Raytracer::applyShading(const Ray& ray, int depth, HitRecord& rec, const SamplingContext& sampling_context, const CameraContext& cam_context) const
{
	Material mat = *rec.material;
	Color color(0.0, 0.0, 0.0);

	if ((mat.type) == "mirror")
	{
		glm::vec3 wo = -ray.direction;
		glm::vec3 wr = (rec.normal * (2 * (glm::dot(rec.normal, wo)))) - wo;
		Ray reflectedRay = Ray(rec.point + rec.normal * render_context.shadow_ray_epsilon, wr, ray.time);

		if (mat.roughness != 0)
		{
			reflectedRay.perturb(mat.roughness);
		}

		color += computeColor(reflectedRay, depth - 1, sampling_context, cam_context) * Color(mat.mirror_reflectance);
	}
	else if ((mat.type) == "conductor")
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

		Ray reflectedRay = Ray(rec.point + rec.normal * render_context.shadow_ray_epsilon, wr, ray.time);
		if (mat.roughness != 0)
		{
			reflectedRay.perturb(mat.roughness);
		}

		color += computeColor(reflectedRay, depth - 1, sampling_context, cam_context) * f_r * mat.mirror_reflectance;
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
		Ray reflectedRay(rec.point + normal * (float)render_context.shadow_ray_epsilon, wr, ray.time);

		reflectedRay.inside = ray.inside;

		if (mat.roughness != 0) reflectedRay.perturb(mat.roughness);

		Color refractedColor(0, 0, 0);
		if (can_refract)
		{
			double cosThetaT = std::sqrt(std::max(0.0, 1.0 - sin2ThetaT));
			glm::vec3 wt = -wo * (float)eta + normal * (float)(eta * cosTheta - cosThetaT);
			wt = glm::normalize(wt);

			Ray refractedRay(rec.point - normal * (float)render_context.shadow_ray_epsilon, wt, ray.time);

			refractedRay.inside = entering;

			if (mat.roughness != 0) refractedRay.perturb(mat.roughness);

			refractedColor = computeColor(refractedRay, depth - 1, sampling_context, cam_context);
		}

		Color reflectedColor = computeColor(reflectedRay, depth - 1, sampling_context, cam_context);
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

	if (!rec.textures.empty())
	{
		//texture_mapping
		for (Texture* tex : rec.textures)
		{
			switch (tex->d_mode)
			{
				case (DecalMode::replace_kd):
				{
					kd = lookupTexture(tex, rec.uv, rec.point);
					break;
				}
				case (DecalMode::blend_kd):
				{
					kd = (lookupTexture(tex, rec.uv, rec.point) + kd) / 2.0f;
					break;
				}
				case (DecalMode::replace_ks):
				{
					ks = lookupTexture(tex, rec.uv, rec.point);
					break;
				}
				case (DecalMode::replace_normal):
				{
					rec.normal = lookupNormalMap(tex, rec);
					break;
				}
				case (DecalMode::bump_normal):
				{
					rec.normal = lookupBumpMap(tex, rec);
					break;
				}
				case (DecalMode::replace_all):
				{
					return lookupTexture(tex, rec.uv, rec.point) * 255.0f;
				}

				default:
				{
					throw std::runtime_error("Invalid decal mode");
				}
			}
		}
	}

	color += Color(ka)* Color(scene->light_sources.ambient_light);

	for (const auto& light : scene->light_sources.point_lights)
	{
		glm::vec3 wi = glm::vec3(light.position) - rec.point;
		double distance = glm::length(wi);
		wi = glm::normalize(wi);
		Ray shadowRay = Ray(rec.point + rec.normal * render_context.shadow_ray_epsilon, wi, ray.time);
		Ray shadowRayPlane = Ray(rec.point
			+ rec.normal
			* static_cast<float>(render_context.shadow_ray_epsilon)
			, wi, ray.time);
		HitRecord shadowRec;
		HitRecord planeShadowRec;
		if (!scene->world->intersect<true>(shadowRay, Interval(render_context.shadow_ray_epsilon, distance), shadowRec)
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
	for (const auto& light : scene->light_sources.area_lights)
	{
		glm::vec3 light_sample_point = (*sampling_context.area_light_samples)[depth - 1][light.id][sampling_context.sample_index];

		// Calculate Shadow Ray towards this specific random point
		glm::vec3 wi = light_sample_point - rec.point;
		double distance = glm::length(wi);
		wi = glm::normalize(wi);
		Ray shadowRay = Ray(rec.point + rec.normal * render_context.shadow_ray_epsilon, wi, ray.time);
		Ray shadowRayPlane = Ray(rec.point
			+ rec.normal
			* static_cast<float>(render_context.shadow_ray_epsilon)
			, wi, ray.time);
		HitRecord shadowRec;
		HitRecord planeShadowRec;
		if (!scene->world->intersect<true>(shadowRay, Interval(render_context.shadow_ray_epsilon, distance), shadowRec)
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

bool Raytracer::hitPlanes(const Ray& ray, Interval ray_t, HitRecord& rec) const
{
	bool hit_anything = false;
	auto closest_so_far = ray_t.max;
	HitRecord temp_rec;

	for (const auto& plane : scene->planes)
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



