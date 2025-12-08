#ifndef SCENE_H
#define SCENE_H

#include <vector>
#include "../camera/pinhole_camera.h"
#include "../camera/distribution_camera.h"
#include "../math_core/math_core.h"
#include "color.h"
#include "../include/parser.hpp"
#include "../material/material_manager.h"
#include "bvh.h"
#include "../light/light.h"
#include "../objects/tlas_box.h"


class Scene{
public:
	Scene();
	Scene(const Scene_& raw_scene, std::vector<TLASBox>& objects);
	~Scene();

	std::vector<PinholeCamera> pinhole_cameras;
	std::vector<DistributionCamera> distribution_cameras;

	Color background_color;
	LightSources light_sources;
	std::shared_ptr<BVH<TLASBox>> world;
};

#endif //SCENE_H
