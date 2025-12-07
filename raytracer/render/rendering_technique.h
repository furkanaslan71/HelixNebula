#ifndef RENDERING_TECHNIQUE_H
#define RENDERING_TECHNIQUE_H

#include "../include/color.h"
#include "../include/ray.h"
#include "../light/light.h"
#include "../material/material_manager.h"
#include "../include/bvh.h"

#define BACKFACE_CULLING 0
typedef struct RendererInfo {
	float shadow_ray_epsilon;
	float intersection_test_epsilon;
	int max_recursion_depth;
	RendererInfo(float s, float i, int d)
	{
		shadow_ray_epsilon = s;
		intersection_test_epsilon = i;
		max_recursion_depth = d;
	}

}RendererInfo;

struct RenderContext {
	const std::vector<std::vector<std::vector<glm::vec3>>>* area_light_samples;

	int sample_index;
};


class RenderingTechnique {
public:
	RenderingTechnique() = default;
	virtual ~RenderingTechnique() = default;
	virtual Color traceRay(const Ray& ray, const RenderContext& context) const = 0;
};

#endif // RENDERING_TECHNIQUE_H

