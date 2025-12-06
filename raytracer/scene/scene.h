#ifndef SCENE_H
#define SCENE_H

#include <vector>
#include "../camera/pinhole_camera.h"
#include "../camera/distribution_camera.h"
#include "vec3.h"
#include "color.h"
#include "../include/parser.hpp"
#include "../material/material_manager.h"
#include "bvh.h"
#include "../light/light.h"


class Scene{
public:
	Scene();
	Scene(const Scene_& raw_scene, std::vector<std::shared_ptr<Hittable>>& objects);
	~Scene();
	glm::mat4 calculateCompositeTransformationMatrix(
		const std::vector<Translation_>& translations,
		const std::vector<Scaling_>& scalings,
		const std::vector<Rotation_>& rotations,
		const PointLight_& raw_light);

	glm::mat4 calculateCompositeTransformationMatrix(
		const std::vector<Translation_>& translations,
		const std::vector<Scaling_>& scalings,
		const std::vector<Rotation_>& rotations,
		const AreaLight_& raw_light);

	std::vector<PinholeCamera> pinhole_cameras;
	std::vector<DistributionCamera> distribution_cameras;

	Color background_color;
	LightSources light_sources;
	std::shared_ptr<BVH> world;
};

#endif //SCENE_H
