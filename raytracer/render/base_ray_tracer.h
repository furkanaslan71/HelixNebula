#ifndef BASE_RAY_TRACER_H
#define BASE_RAY_TRACER_H
#include "rendering_technique.h"
#include "../objects/plane.h"
#include "../random/my_random.h"
#include "../objects/tlas_box.h"
#define OUT

class BaseRayTracer : public RenderingTechnique {
public:
	BaseRayTracer( Color& background_color,
		LightSources& light_sources,
		std::shared_ptr<BVH<TLASBox>> world,
		std::vector<Plane>& planes,
		MaterialManager& material_manager,
		RendererInfo& renderer_info);

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
};

#endif // BASE_RAY_TRACER_H