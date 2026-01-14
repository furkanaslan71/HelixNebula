#include "raytracer.h"
#include <typeinfo>

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

Raytracer::Raytracer(std::unique_ptr<Scene> _scene, const RenderContext& _render_context)
    : render_context(_render_context)
{
    scene = std::move(_scene);
}

void Raytracer::renderScene()
{
    std::string saveDir = FS::absolute(__FILE__).parent_path() / "../../outputs";
    for (const auto& cam : scene->cameras)
    {
        std::vector<std::vector<Color>> image;
        renderOneCamera(cam, image);
    	std::string image_name = cam->image_name;
    	std::string extension = getFileExtension(image_name);

    	if (cam->flip_x)
    	{
    		for (auto& row : image)
    		{
    			std::reverse(row.begin(), row.end());
    		}
    	}

    	if (extension== "png" || extension == "jpg" || extension == ".jpeg")
			saveImage(saveDir, image_name, image, ImageType::SDR);
    	else if (extension == "exr" || extension == "hdr")
    	{

    		// Check if we should save the raw EXR based on filename suffix
    		bool shouldSaveRaw = false;
    		if (image_name.length() >= 4) {
    			std::string lower_name = image_name;
    			// Simple suffix check
    			if (lower_name.ends_with(".exr") || lower_name.ends_with(".hdr") ||
					lower_name.ends_with("_exr") || lower_name.ends_with("_hdr")) {
    				shouldSaveRaw = true;
					}
    		}

    		if (shouldSaveRaw) {
    			std::string raw_path = image_name;
    			// Ensure it ends with .exr for the filesystem
    			if (!raw_path.ends_with(".exr")) {
    				raw_path += ".exr";
    			}
    			saveEXR(saveDir + "/" + raw_path, image);
    		}
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

void Raytracer::renderOneCamera(std::shared_ptr<BaseCamera> camera, std::vector<std::vector<Color>>& output)
{
	render_context.max_recursion_depth = camera->recursion_depth;
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
			if (render_context.env_map != nullptr)
			{
				return Color(lookupEnvMap(scene->light_sources.env_light,
				   glm::normalize(ray.direction)));
			}

			if (render_context.b_type == BackgroundType::Color)
				return Color(render_context.background_info.background_color);

			return Color(lookupBackgroundTex(render_context.background_info.background_tex,
			   glm::normalize(ray.direction), cam_context));
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

		Color reflectedColor = computeColor(reflectedRay, depth - 1, sampling_context, cam_context);

		color += reflectedColor * Color(mat.mirror_reflectance);
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
					auto res = Color(lookupTexture(tex, rec.uv, rec.point, true));
					return  res;
				}

				default:
				{
					throw std::runtime_error("Invalid decal mode");
				}
			}
		}
	}

	if (mat.degamma)
	{
		kd = glm::pow(kd, glm::vec3(2.2));
		ks = glm::pow(ks, glm::vec3(2.2));
		ka = glm::pow(ka, glm::vec3(2.2));
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
	for (const auto& light : scene->light_sources.directional_lights)
	{
		glm::vec3 wi = -glm::normalize(light.direction);
		Ray shadowRay = Ray(rec.point + rec.normal * render_context.shadow_ray_epsilon, wi, ray.time);
		Ray shadowRayPlane = Ray(rec.point
			+ rec.normal
			* static_cast<float>(render_context.shadow_ray_epsilon)
			, wi, ray.time);
		HitRecord shadowRec;
		HitRecord planeShadowRec;
		if (!scene->world->intersect<true>(shadowRay, Interval(render_context.shadow_ray_epsilon, INFINITY), shadowRec)
			&& !hitPlanes(shadowRayPlane, Interval(0, INFINITY), planeShadowRec))
		{
			//diffuse
			color += Color(kd * light.radiance * glm::max(0.0f, glm::dot(rec.normal, wi)));
			//specular
			glm::vec3 wo = (ray.origin - rec.point);
			wo = glm::normalize(wo);

			glm::vec3 h = glm::normalize(wi + wo);
			float specular = std::pow(glm::max(0.0f, glm::dot(rec.normal, h)),
									  mat.phong_exponent);
			color += Color(ks * light.radiance * specular);
		}
	}
	for (const auto& light : scene->light_sources.spot_lights)
	{
	    glm::vec3 wi = light.position - rec.point;
	    double distance = glm::length(wi);
	    wi = glm::normalize(wi);

	    Ray shadowRay(
	        rec.point + rec.normal * render_context.shadow_ray_epsilon,
	        wi,
	        ray.time
	    );

	    Ray shadowRayPlane(
	        rec.point + rec.normal * render_context.shadow_ray_epsilon,
	        wi,
	        ray.time
	    );

	    HitRecord shadowRec;
	    HitRecord planeShadowRec;

	    if (!scene->world->intersect<true>(
	            shadowRay,
	            Interval(render_context.shadow_ray_epsilon, distance),
	            shadowRec)
	        && !hitPlanes(shadowRayPlane, Interval(0, distance), planeShadowRec))
	    {
	        glm::vec3 spotDir = glm::normalize(light.direction);
	        double cosAlpha = glm::dot(-wi, spotDir);

	        double cosCoverage = std::cos(glm::radians(light.coverage_angle * 0.5f));
	        double cosFalloff  = std::cos(glm::radians(light.falloff_angle  * 0.5f));

	        double spotFactor = 0.0;

	        if (cosAlpha >= cosFalloff)
	        {
	            spotFactor = 1.0;
	        }
	        else if (cosAlpha >= cosCoverage)
	        {
	            double s = (cosAlpha - cosCoverage) /
	                       (cosFalloff - cosCoverage);
	            spotFactor = std::pow(s, 4.0);
	        }
	        else
	        {
	            continue;
	        }

	        double attenuation = spotFactor / (distance * distance);

	        double NdotL = std::max(0.0f, glm::dot(rec.normal, wi));
	        color += Color(kd) * Color(light.intensity) * attenuation * NdotL;

	        glm::vec3 wo = glm::normalize(ray.origin - rec.point);
	        glm::vec3 h  = glm::normalize(wi + wo);

	        double spec = std::pow(
	            std::max(0.0f, glm::dot(rec.normal, h)),
	            mat.phong_exponent
	        );

	        color += Color(ks) * Color(light.intensity) * attenuation * spec;
	    }
	}

	if (scene->light_sources.env_light.img != nullptr && mat.type != "mirror")
	{
		glm::vec3 n = rec.normal;
		glm::vec3 u = rec.surface_tangents.u;
		glm::vec3 v = rec.surface_tangents.v;
		float pi = glm::pi<float>();

		Color env_contribution(0, 0, 0);
		int samples = 25;

		for (int s = 0; s < samples; ++s)
		{
			float s1 = generateRandomFloat(0, 1);
			float s2 = generateRandomFloat(0, 1);
			glm::vec3 l;
			Color sample_color(0, 0, 0);

			if (scene->light_sources.env_light.sampler == Sampler::uniform)
			{
				l = u * glm::sqrt(1 - s1 * s1) * glm::cos(s2 * 2 * pi) +
					v * glm::sqrt(1 - s1 * s1) * glm::sin(s2 * 2 * pi) + n * s1;

				l = glm::normalize(l);

				float cosTheta = std::max(0.0f, glm::dot(n, l));
				sample_color = Color(lookupEnvMap(scene->light_sources.env_light, l)) * (2.0f * pi * cosTheta);
			}
			else if (scene->light_sources.env_light.sampler == Sampler::cosine)
			{
				l = u * glm::sqrt(s1) * glm::cos(s2 * 2 * pi) +
					v * glm::sqrt(s1) * glm::sin(s2 * 2 * pi) + n * glm::sqrt(1 - s1);

				l = glm::normalize(l);

				float theta = glm::asin(glm::sqrt(s1));
				sample_color = Color(lookupEnvMap(scene->light_sources.env_light, l) * pi / glm::cos(theta));
			}

			Ray envShadowRay(rec.point + n * (float)render_context.shadow_ray_epsilon, l, ray.time);
			HitRecord envRec;

			if (!scene->world->intersect<true>(envShadowRay, Interval(render_context.shadow_ray_epsilon, INFINITY), envRec)
				&& !hitPlanes(envShadowRay, Interval(render_context.shadow_ray_epsilon, INFINITY), envRec))
			{
				env_contribution += sample_color;
			}
		}

		color += env_contribution / (float)samples;
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



