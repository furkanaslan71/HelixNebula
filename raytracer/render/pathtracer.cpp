//
// Created by furkan on 1/14/26.
//
#include "pathtracer.h"

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


PathTracer::PathTracer(std::unique_ptr<Scene> _scene, const RenderContext& _render_context)
    : render_context(_render_context)
{
    scene = std::move(_scene);
}

Color PathTracer::tracePath(const Ray &initial_ray, const CameraContext &cam_context) const
{
    Color accumulated_light(0.0f);
    Color throughput(1.0f);
    Ray ray = initial_ray;

    for (int depth = 0; depth <= render_context.max_recursion_depth; ++depth)
    {
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
                    accumulated_light += throughput * lookupEnvMap(scene->light_sources.env_light,
                       glm::normalize(ray.direction));
                }

                else if (render_context.b_type == BackgroundType::Color)
                {
                    accumulated_light += throughput * render_context.background_info.background_color;
                }
                else
                {
                    accumulated_light += throughput * lookupBackgroundTex(render_context.background_info.background_tex,
                   glm::normalize(ray.direction), cam_context);
                }
                break;
            }
        }

        if (rec.radiance.has_value())
        {
            if (glm::dot(rec.normal, -ray.direction) > 0)
            {
                accumulated_light += throughput * Color(rec.radiance.value());
            }
            else
            {
                ray = Ray(rec.point + ray.direction * (float)render_context.shadow_ray_epsilon, ray.direction, ray.time);
                continue;
            }
        }

        Material mat = *rec.material;
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

        Ray next_ray;
        if ((mat.type) == "mirror")
	    {
		    glm::vec3 wo = -ray.direction;
		    glm::vec3 wr = (rec.normal * (2 * (glm::dot(rec.normal, wo)))) - wo;
		    next_ray = Ray(rec.point + rec.normal * render_context.shadow_ray_epsilon, wr, ray.time);

		    if (mat.roughness != 0)
		    {
			    next_ray.perturb(mat.roughness);
		    }

		    throughput *= Color(mat.mirror_reflectance);
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

		    next_ray = Ray(rec.point + rec.normal * render_context.shadow_ray_epsilon, wr, ray.time);
		    if (mat.roughness != 0)
		    {
			    next_ray.perturb(mat.roughness);
		    }

		    throughput *= Color((float)f_r * mat.mirror_reflectance);
	    }
	    else if (mat.type == "dielectric")
        {
            glm::vec3 wo = -ray.direction;
            glm::vec3 geometric_normal = rec.normal;
            bool entering = glm::dot(ray.direction, geometric_normal) < 0;

            double n1, n2;
            glm::vec3 normal;
            double current_ior = (double)mat.refraction_index;

            if (entering) { n1 = 1.0; n2 = current_ior; normal = geometric_normal; }
            else { n1 = current_ior; n2 = 1.0; normal = -geometric_normal; }

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

            // Monte Carlo Path Selection: Reflect or Refract?
            float rnd = generateRandomFloat(0, 1);

            if (rnd < F_r)
            {
                // -- REFLECTION --
                glm::vec3 wr = reflect(wo, normal);
                next_ray = Ray(rec.point + normal * (float)render_context.shadow_ray_epsilon, wr, ray.time);

                // Maintain previous inside state
                next_ray.inside = ray.inside;

                if (mat.roughness != 0) next_ray.perturb(mat.roughness);

                // Throughput logic: multiply by F_r / P_reflect.
                // Since P_reflect = F_r, they cancel out. throughput *= 1.0;
            }
            else
            {
                // -- REFRACTION --
                double cosThetaT = std::sqrt(std::max(0.0, 1.0 - sin2ThetaT));
                glm::vec3 wt = -wo * (float)eta + normal * (float)(eta * cosTheta - cosThetaT);
                wt = glm::normalize(wt);

                next_ray = Ray(rec.point - normal * (float)render_context.shadow_ray_epsilon, wt, ray.time);

                // Update inside state
                next_ray.inside = entering;

                if (mat.roughness != 0) next_ray.perturb(mat.roughness);

                // Beer's Law Absorption (if ray was traveling inside the medium)
                if (!entering)
                {
                     double dist = rec.t;
                     glm::vec3 absorb;
                     absorb.x = std::exp(-mat.absorption_coefficient.x * dist);
                     absorb.y = std::exp(-mat.absorption_coefficient.y * dist);
                     absorb.z = std::exp(-mat.absorption_coefficient.z * dist);
                     throughput *= Color(absorb);
                }

                // Throughput logic: multiply by (1-F_r) / P_refract.
                // Since P_refract = (1-F_r), they cancel out. throughput *= 1.0;
            }
        }
	    else
	    {
	        // 1. Build Local Coordinate System (ONB)
	        glm::vec3 w = glm::normalize(rec.normal);
	        glm::vec3 a = (std::abs(w.x) > 0.9f) ? glm::vec3(0, 1, 0) : glm::vec3(1, 0, 0);
	        glm::vec3 v = glm::normalize(glm::cross(w, a));
	        glm::vec3 u = glm::cross(v, w);

	        // 2. Cosine Weighted Hemisphere Sampling (Samples wi)
	        float r1 = generateRandomFloat(0, 1);
	        float r2 = generateRandomFloat(0, 1);
	        float r = std::sqrt(r1);
	        float theta = 2.0f * M_PI * r2;

	        float x = r * std::cos(theta);
	        float y = r * std::sin(theta);
	        float z = std::sqrt(std::max(0.0f, 1.0f - r1)); // z is cosTheta

	        glm::vec3 local_dir(x, y, z);
	        glm::vec3 next_dir = glm::normalize(x * u + y * v + z * w);

	        next_ray = Ray(rec.point + next_dir * (float)render_context.shadow_ray_epsilon, next_dir, ray.time);

	        // 3. Evaluate BRDF
	        // wo = -ray.direction (Unit vector towards eye)
	        // wi = next_dir (Unit vector towards light/next bounce)
	        glm::vec3 wo = -glm::normalize(ray.direction);

	        Color brdf_val(0.0f);

	        // Check if material has a specific BRDF definition
	        if (rec.material->brdf) {
	            // Pass kd and ks from material (which might have been texturized above)
	            brdf_val = Color(rec.material->brdf->Evaluate(next_dir, wo, rec.normal, kd, ks));
	        }
	        else {
	            // Fallback to Lambertian Diffuse if no specific BRDF is set
	            // f = kd / PI
	            brdf_val = Color(kd * (float)(1.0 / M_PI));
	        }

	        // 4. Update Throughput
	        // Formula: throughput * (BRDF * cosTheta) / PDF
	        // PDF (Cosine Sampling) = cosTheta / PI
	        // Factor = BRDF * PI

	        throughput *= brdf_val * (float)M_PI;
	    }
        ray = next_ray;

        // Terminate if throughput is black to save performance
        if (throughput.r <= 1e-5 && throughput.g <= 1e-5 && throughput.b <= 1e-5) break;
    }
    return accumulated_light;
}


void PathTracer::renderScene()
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

void PathTracer::renderOneCamera(std::shared_ptr<BaseCamera> camera, std::vector<std::vector<Color>>& output)
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
        threads[threadId] = std::thread(&PathTracer::renderLoop, this, threadId, numThreads, camera, std::ref(output));
    }
    for (auto& t : threads) t.join();
#endif
}

void PathTracer::generateRaySamples(std::shared_ptr<BaseCamera> camera, int i, int j,
                        std::vector<std::pair<glm::vec3, glm::vec3>>& out_samples,
                        std::mt19937& rng) const
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

void PathTracer::renderLoop(int threadId, int stride, std::shared_ptr<BaseCamera> camera, std::vector<std::vector<Color>>& output) const
{
    std::vector<std::pair<glm::vec3, glm::vec3>> ray_samples;
    ray_samples.reserve(camera->num_samples);


    std::mt19937 rng(std::random_device{}());

    for (int i = threadId; i < camera->image_height; i += stride)
    {
        for (int j = 0; j < camera->image_width; ++j)
        {
#if DEBUG_PIXEL
            if (i != HEIGHT || j != WIDTH) continue; // Only process the debug pixel
#endif

            generateRaySamples(camera, i, j, ray_samples, rng);

            Color pixel_color(0.0, 0.0, 0.0);

            for (int k = 0; k < camera->num_samples; k++)
            {
                Ray primary_ray(ray_samples[k].first, ray_samples[k].second, generateRandomFloat(0, 1));
                pixel_color += tracePath(primary_ray, camera->context);
            }
            output[i][j] = (pixel_color / (float)camera->num_samples);
        }
    }
}


bool PathTracer::hitPlanes(const Ray& ray, Interval ray_t, HitRecord& rec) const
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