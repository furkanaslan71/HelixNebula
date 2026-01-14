#ifndef RAYTRACER_H
#define RAYTRACER_H

#include <vector>
#include "config.h"
#include "core/hittable.h"
#include "scene/scene.h"
#include <filesystem>
#include "io/image_io.h"
#include "texture_mapping/texture_lookup.h"
#include "rendercontext.h"



struct SamplingContext {
    const std::vector<std::vector<std::vector<glm::vec3>>>* area_light_samples;
    int sample_index;
};

class Raytracer {
public:
    Raytracer(std::unique_ptr<Scene> _scene, const RenderContext& _render_context);

    void renderScene();
private:
    void renderOneCamera(std::shared_ptr<BaseCamera> camera, std::vector<std::vector<Color>>& output);

    Color traceRay(const Ray& ray, const SamplingContext& sampling_context, const CameraContext& cam_context) const ;

    Color computeColor(const Ray& ray, int depth, const SamplingContext& sampling_context, const CameraContext& cam_context) const;

    Color applyShading(const Ray& ray, int depth, HitRecord& rec, const SamplingContext& sampling_context, const CameraContext& cam_context) const;

    bool hitPlanes(const Ray& ray, Interval ray_t, HitRecord& rec) const;

    void renderLoop(int threadId, int stride, std::shared_ptr<BaseCamera> camera, std::vector<std::vector<Color>>& output) const;

    std::unique_ptr<Scene> scene;
    RenderContext render_context;
};


#endif //RAYTRACER_H
