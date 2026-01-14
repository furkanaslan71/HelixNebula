//
// Created by furkan on 1/14/26.
//

#ifndef PATHTRACER_H
#define PATHTRACER_H
#include "rendercontext.h"
#include "scene/scene.h"
#include "camera/base_camera.h"

class PathTracer {
public:
    PathTracer(std::unique_ptr<Scene> _scene, const RenderContext& _render_context);

    void renderScene();
private:
    void renderOneCamera(std::shared_ptr<BaseCamera> camera, std::vector<std::vector<Color>>& output);

    void generateRaySamples(std::shared_ptr<BaseCamera> camera, int i, int j,
                                   std::vector<std::pair<glm::vec3, glm::vec3>>& out_samples,
                                   std::mt19937& rng) const;

    Color tracePath(const Ray& ray, const CameraContext& cam_context) const;

    bool hitPlanes(const Ray& ray, Interval ray_t, HitRecord& rec) const;

    void renderLoop(int threadId, int stride, std::shared_ptr<BaseCamera> camera, std::vector<std::vector<Color>>& output) const;

    std::unique_ptr<Scene> scene;
    RenderContext render_context;
};

#endif //PATHTRACER_H
