#ifndef BASE_CAMERA_H
#define BASE_CAMERA_H
#include <vector>
#include <thread>

#include "config.h"
#include "core/color.h"
#include "core/ray.h"
#include "math_core/math_core.h"
#include "acceleration/bvh.h"
#include "parser/parser.hpp"
#include "light/light.h"


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

	std::string image_name;
	int num_samples;
	int id;


	void generatePixelSamples(int i, int j, std::vector<glm::vec3>& out_samples) const;

	std::vector<AreaLight>& area_lights;
	
	glm::vec3 w, u, v, q, su, sv, m;
	glm::vec3 position, gaze, up;
	double near_plane[4]; // l, r, b, t
	double near_distance;
	int image_width, image_height;
	int recursion_depth;
	int num_area_lights;
	CameraContext context;

};

#endif //BASE_CAMERA_H
