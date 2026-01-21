//
// Created by furkan on 1/14/26.
//
#include "pathtracer.h"

glm::vec3 isotropicSample(float r1, float r2)
{
    float phi = 2.0f * M_PI * r1;
    float cosTheta = 1.0f - 2.0f * r2;
    float sinTheta = std::sqrt(std::max(0.0f, 1.0f - cosTheta * cosTheta));
    return glm::vec3(sinTheta * std::cos(phi), sinTheta * std::sin(phi), cosTheta);
}

glm::vec3 sampleHenyeyGreenstein(glm::vec3 incoming_dir, float g, float r1, float r2)
{
    // If g is very small, treat as isotropic to avoid division by zero
    if (std::abs(g) < 1e-3)
    {
        return isotropicSample(r1, r2);
    }

    // 1. Sample cosTheta using the HG formula
    float sqrTerm = (1.0f - g * g) / (1.0f - g + 2.0f * g * r1);
    float cosTheta = (1.0f + g * g - sqrTerm * sqrTerm) / (2.0f * g);
    float sinTheta = std::sqrt(std::max(0.0f, 1.0f - cosTheta * cosTheta));
    float phi = 2.0f * M_PI * r2;

    // 2. Create local direction (z-axis is the "forward" incoming direction)
    glm::vec3 local_dir(sinTheta * std::cos(phi), sinTheta * std::sin(phi), cosTheta);

    // 3. Align local_dir to the incoming_dir (World Space transformation)
    glm::vec3 w = incoming_dir; // This is our 'z'
    glm::vec3 a = (std::abs(w.x) > 0.9f) ? glm::vec3(0, 1, 0) : glm::vec3(1, 0, 0);
    glm::vec3 v = glm::normalize(glm::cross(w, a));
    glm::vec3 u = glm::cross(v, w);

    return glm::normalize(local_dir.x * u + local_dir.y * v + local_dir.z * w);
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

PathTracer::PathTracer(std::unique_ptr<Scene> _scene, const RenderContext& _render_context)
    : render_context(_render_context)
{
    scene = std::move(_scene);
}

Color PathTracer::tracePath(const Ray &initial_ray, const CameraContext &cam_context) const
{
    Color accumulated_light(0.0f);
    int max_rec_depth = render_context.max_recursion_depth;
    std::vector<PathState> stack;
    stack.reserve(cam_context.splitting_factor * 2 + max_rec_depth);
    stack.push_back({initial_ray, Color(1.0f), 0});

    bool russian_roulette = cam_context.options.at(RendererParams::RussianRoulette);
    float survival_rate;
    if (russian_roulette)  max_rec_depth = 16;

    while (!stack.empty())
    {
        PathState state = stack.back();
        stack.pop_back();

        // Check recursion limit
        if (state.depth > max_rec_depth) continue;

        Ray ray = state.ray;
        Color throughput = state.throughput;

        if (russian_roulette && state.depth >= cam_context.min_recursion_depth)
        {
            survival_rate = std::clamp(throughput.max(), 1e-7, 0.99);

            float bullet_to_the_head = generateRandomFloat(0, 1);
            if (bullet_to_the_head > survival_rate)
                continue;

            throughput *= 1.0f / survival_rate;
        }

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
                continue;
            }
        }

        if (ray.inside)
        {
            Material mat = *rec.material;
            glm::vec3 sigma_s = mat.scattering_coefficient;
            glm::vec3 sigma_a = mat.absorption_coefficient;
            glm::vec3 sigma_t = sigma_a + sigma_s;
            //float max_sigma_t = std::max({sigma_t.x, sigma_t.y, sigma_t.z});
            auto luminance = [](const glm::vec3& a)
            {
                return 0.2126 * a.r + 0.7152 * a.g + 0.0722 * a.b;
            };
            float lum_sigma_t = std::max(1e-6, luminance(sigma_t));

            if (lum_sigma_t < 0)
                throw std::runtime_error("no absorption and scattering values");

            float distance = -std::log(1 - generateRandomFloat(0, 1)) / lum_sigma_t;
            if (distance > rec.t)
            {
                float d = rec.t;
                glm::vec3 Tr = glm::exp(-sigma_t * d);
                //float pdf = lum_sigma_t * std::exp(-lum_sigma_t * d);
                float pdf = std::exp(-lum_sigma_t * d);
                //pdf = 1.0;
                throughput *= Tr / pdf;
                ray.origin += d * ray.direction;
                ray.inside = false;
                rec.material->type = "none";
            }
            else
            {
                ray.origin += distance * ray.direction;
                float pdf = lum_sigma_t * std::exp(-lum_sigma_t * distance);

                throughput *= sigma_s * glm::exp(-sigma_t * distance) / pdf;

                float scatter_prob = luminance(sigma_s) / lum_sigma_t;
                if (generateRandomFloat(0, 1) > scatter_prob)
                    continue;

                ray.direction = glm::normalize(isotropicSample(generateRandomFloat(0, 1),
                    generateRandomFloat(0, 1)));

                throughput *= 1.0f / scatter_prob;

                stack.push_back({ray, throughput, state.depth + 1});
                continue;
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
                Ray pass_ray(rec.point + ray.direction * (float)render_context.shadow_ray_epsilon, ray.direction, ray.time);
                stack.push_back({pass_ray, throughput, state.depth}); // Same depth, effectively skipping the object
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
                        return Color(lookupTexture(tex, rec.uv, rec.point, true));
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

        int num_splits = 1;
        if (state.depth == 0 && cam_context.splitting_factor > 1)
        {
            num_splits = cam_context.splitting_factor;
        }


        for (int i = 0; i < num_splits; ++i)
        {
            Color next_throughput = throughput;
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

                next_throughput *= Color(mat.mirror_reflectance);
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

                next_throughput *= Color((float)f_r * mat.mirror_reflectance);
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
                    //std::cout << F_r << std::endl;
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
                    if (false && !entering)
                    {
                        double dist = rec.t;
                        glm::vec3 absorb;
                        absorb.x = std::exp(-mat.absorption_coefficient.x * dist);
                        absorb.y = std::exp(-mat.absorption_coefficient.y * dist);
                        absorb.z = std::exp(-mat.absorption_coefficient.z * dist);
                        next_throughput *= Color(absorb);
                    }

                    // Throughput logic: multiply by (1-F_r) / P_refract.
                    // Since P_refract = (1-F_r), they cancel out. throughput *= 1.0;
                }
            }
            else
            {
                glm::vec3 w = glm::normalize(rec.normal);
                glm::vec3 a = (std::abs(w.x) > 0.9f) ? glm::vec3(0, 1, 0) : glm::vec3(1, 0, 0);
                glm::vec3 v = glm::normalize(glm::cross(w, a));
                glm::vec3 u = glm::cross(v, w);

                glm::vec3 next_dir;
                float r1 = generateRandomFloat(0, 1);
                float r2 = generateRandomFloat(0, 1);
                float theta = 2.0f * M_PI * r2;
                float r, x, y, z;
                float invPDF_cosTheta;

                if (cam_context.options.at(RendererParams::ImportanceSampling))
                {
                    r = std::sqrt(r1);
                    x = r * std::cos(theta);
                    y = r * std::sin(theta);
                    z = std::sqrt(std::max(0.0f, 1.0f - r1)); // z is cosTheta
                    next_dir = glm::normalize(x * u + y * v + z * w);

                    invPDF_cosTheta = M_PI;
                }
                else
                {
                    z = r1;
                    r = std::sqrt(std::max(0.0f, 1.0f - z * z));
                    x = r * std::cos(theta);
                    y = r * std::sin(theta);
                    next_dir = glm::normalize(x * u + y * v + z * w);

                    float cosTheta = std::max(0.0f, glm::dot(rec.normal, next_dir));
                    invPDF_cosTheta = cosTheta * 2.0f * M_PI;
                }

                next_ray = Ray(rec.point + next_dir * (float)render_context.shadow_ray_epsilon, next_dir, ray.time);
                glm::vec3 wo = -glm::normalize(ray.direction);

                Color brdf_val(0.0f);
                if (rec.material->brdf)
                {
                    brdf_val = Color(rec.material->brdf->Evaluate(next_dir, wo, rec.normal, kd, ks));
                }
                else
                {
                    // Fallback to Lambertian Diffuse if no specific BRDF is set
                    // f = kd / PI
                    brdf_val = Color(kd * (float)(1.0 / M_PI));
                }
                next_throughput *= brdf_val * invPDF_cosTheta;
            }

            if (num_splits > 1)
            {
                next_throughput = next_throughput / num_splits;
            }

            if (next_throughput.r > 1e-5 || next_throughput.g > 1e-5 || next_throughput.b > 1e-5)
            {
                stack.push_back({next_ray, next_throughput, state.depth + 1});
            }
        }
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
        //break;
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