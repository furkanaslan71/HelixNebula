#ifndef SCENE_H
#define SCENE_H
#include <vector>

#include "camera/pinhole_camera.h"
#include "camera/distribution_camera.h"
#include "math_core/math_core.h"
#include "core/color.h"
#include "parser/parser.hpp"
#include "acceleration/bvh.h"
#include "light/light.h"
#include "objects/tlas_box.h"
#include "objects/plane.h"

class Scene{
public:
	Scene();
	Scene(const Scene_& raw_scene, std::vector<TLASBox>& objects, const std::vector<Plane>& _planes);
	~Scene();

	//std::vector<PinholeCamera> pinhole_cameras;
	//std::vector<DistributionCamera> distribution_cameras;
	std::vector<std::shared_ptr<BaseCamera>> cameras;

	Color background_color;
	LightSources light_sources;
	std::shared_ptr<BVH<TLASBox>> world;
	std::vector<Plane> planes;
};

#endif //SCENE_H
