#ifndef RENDER_MANAGER_H
#define RENDER_MANAGER_H

#include "core/color.h"
#include "scene/scene.h"
#include "base_ray_tracer.h"
#include "material/material_manager.h"
#include "parser/parser.hpp"

class RenderManager {
public:
  RenderManager(const Scene& _scene,
    const MaterialManager& _material_manager,
    const RendererInfo _renderer_info,
    BaseRayTracer& _rendering_technique);
    void render() const;

private:
  void saveImage(const std::string& outputDir, const std::string& fileName, std::vector<std::vector<Color>>& image)
    const;

  const Scene& scene;
  const MaterialManager& material_manager;
  const RendererInfo renderer_info;
	BaseRayTracer& technique;
  
};

#endif // RENDER_MANAGER_H