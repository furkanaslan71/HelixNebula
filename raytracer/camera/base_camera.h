#ifndef BASE_CAMERA_H
#define BASE_CAMERA_H

#define OUT
#define IN


#include <vector>
#include "color.h"
#include "ray.h"
#include "../math_core/math_core.h"
#include "bvh.h"
#include "parser.hpp"
#include "../render/rendering_technique.h"
#include <thread>
#include "../render/base_ray_tracer.h"
#include "../light/light.h"


class Scene;
class RenderingTechnique;

class BaseCamera {
public:
	BaseCamera(
		const Camera_& cam,
		const std::vector<Translation_>& translations,
		const std::vector<Scaling_>& scalings,
		const std::vector<Rotation_>& rotations,
		int _recursion_depth,
		int _num_area_lights,
		std::vector<AreaLight>& _area_lights
	);
	virtual ~BaseCamera() = default;
	virtual void render(IN const BaseRayTracer& rendering_technique,
		OUT std::vector<std::vector<Color>>& image) const = 0;

	std::string image_name;
	int num_samples;
	int id;
	
protected:
	void generatePixelSamples(int i, int j, std::vector<glm::vec3>& out_samples) const;

	std::vector<AreaLight>& area_lights;
	glm::vec3 w, u, v, q, su, sv, m;
	glm::vec3 position, gaze, up;
	double near_plane[4]; // l, r, b, t
	double near_distance;
	int image_width, image_height;
	int recursion_depth;
	int num_area_lights;
};

#endif //BASE_CAMERA_H
