//
// Created by furkan on 1/14/26.
//
#include "pathtracer.h"


PathTracer::PathTracer(std::unique_ptr<Scene> _scene, const RenderContext& _render_context)
    : render_context(_render_context)
{
    scene = std::move(_scene);
}

Color PathTracer::tracePath(const Ray &ray, const CameraContext &cam_context) const
{
    return 0;
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