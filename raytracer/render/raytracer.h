#ifndef RAYTRACER_H
#define RAYTRACER_H

#include <vector>
#include "config.h"
#include "core/hittable.h"
#include "scene/scene.h"
#include <filesystem>
#include "io/save_image.h"


struct RenderContext {
    Color background_color;
    float shadow_ray_epsilon;
    float intersection_test_epsilon;
    int max_recursion_depth;
    RenderContext(Color background_color, float s, float i, int d)
    {
        this->background_color = background_color;
        shadow_ray_epsilon = s;
        intersection_test_epsilon = i;
        max_recursion_depth = d;
    }
};

struct SamplingContext {
    const std::vector<std::vector<std::vector<glm::vec3>>>* area_light_samples;
    int sample_index;
};

class Raytracer {
public:
    Raytracer(std::unique_ptr<Scene> _scene, const RenderContext& _render_context);

    void renderScene() const;
private:
    void renderOneCamera(std::shared_ptr<BaseCamera> camera, std::vector<std::vector<Color>>& output) const;

    Color traceRay(const Ray& ray, const SamplingContext& sampling_context) const ;

    Color computeColor(const Ray& ray, int depth, const SamplingContext& sampling_context) const;

    Color applyShading(const Ray& ray, int depth, HitRecord& rec, const SamplingContext& sampling_context) const;

    bool hitPlanes(const Ray& ray, Interval ray_t, HitRecord& rec) const;

    void renderLoop(int threadId, int stride, std::shared_ptr<BaseCamera> camera, std::vector<std::vector<Color>>& output) const;

    std::unique_ptr<Scene> scene;
    RenderContext render_context;
};


#endif //RAYTRACER_H
