#ifndef BASE_RAY_TRACER_H
#define BASE_RAY_TRACER_H

#include "glm_config.h"
#include "rendering_technique.h"
#include "objects/plane.h"
#include "math_core/math_core.h"
#include "objects/tlas_box.h"
#include "texture_mapping/texture_fetcher.h"


class BaseRayTracer : public RenderingTechnique {
public:
	BaseRayTracer( Color& background_color,
		LightSources& light_sources,
		std::shared_ptr<BVH<TLASBox>> world,
		std::vector<Plane>& planes,
		MaterialManager& material_manager,
		RendererInfo& renderer_info,
		TextureFetcher& texture_fetcher);

	Color traceRay(const Ray& ray, const RenderContext& context) const override;

	Color computeColor(const Ray& ray, int depth, const RenderContext& context) const;

	Color applyShading(const Ray& ray, int depth, HitRecord& rec, const RenderContext& context) const;

	bool hitPlanes(const Ray& ray, Interval ray_t, HitRecord& rec) const;

	Color& background_color;
	LightSources& light_sources;
	std::shared_ptr<BVH<TLASBox>> world;
	std::vector<Plane>& planes;
	MaterialManager& material_manager;
	RendererInfo& renderer_info;
	TextureFetcher& texture_fetcher;
	BackgroundCameraData bgc_data;
	
};

#endif // BASE_RAY_TRACER_H