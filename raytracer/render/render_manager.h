#ifndef RENDER_MANAGER_H
#define RENDER_MANAGER_H

#include "core/color.h"
#include "scene/scene.h"
#include "base_ray_tracer.h"
#include "parser/parser.hpp"

class RenderManager {
public:
  RenderManager(const Scene& _scene,
    const RendererInfo _renderer_info,
    BaseRayTracer& _rendering_technique);
    void render() const;

private:
  void saveImage(const std::string& outputDir, const std::string& fileName, std::vector<std::vector<Color>>& image)
    const;

  const Scene& scene;
  const RendererInfo renderer_info;
	BaseRayTracer& technique;
  
};

#endif // RENDER_MANAGER_H