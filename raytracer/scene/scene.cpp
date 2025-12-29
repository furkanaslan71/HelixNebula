#include "scene.h"


Scene::Scene() {}

Scene::Scene(const Scene_& raw_scene, std::vector<TLASBox>& objects, const std::vector<Plane>& _planes)
	: background_color(raw_scene.background_color.x, raw_scene.background_color.y, raw_scene.background_color.z)
{
	for (const auto& raw_light : raw_scene.point_lights) 
	{
		light_sources.point_lights.push_back(PointLight(raw_light.id,
			glm::vec3(raw_light.position),
			Color(raw_light.intensity.x, raw_light.intensity.y, raw_light.intensity.z)));
	}
	for (const auto& raw_light : raw_scene.area_lights)
	{
		glm::vec3 u, v;

		createONB(raw_light.normal, u, v);

		light_sources.area_lights.push_back(
			AreaLight(
				raw_light.id,
			glm::vec3(raw_light.position),
			glm::vec3(raw_light.normal),
				raw_light.edge,
			glm::vec3(raw_light.radiance.x, raw_light.radiance.y, raw_light.radiance.z),
				u,
				v
			)
		);
	}
	
	light_sources.ambient_light = Color(raw_scene.ambient_light.x,
		raw_scene.ambient_light.y,
		raw_scene.ambient_light.z);

	int recursion_depth = raw_scene.max_recursion_depth;
	int num_area_lights = raw_scene.area_lights.size();

	for (const auto& raw_camera : raw_scene.cameras)
	{
		if (raw_camera.aperture_size == 0)
		{
			cameras.emplace_back(std::make_shared<PinholeCamera>(raw_camera,
				raw_scene.translations,
				raw_scene.scalings,
				raw_scene.rotations,
				recursion_depth,
				num_area_lights,
				this->light_sources.area_lights
				));
		}
		else
		{
			cameras.emplace_back(std::make_shared<DistributionCamera>(raw_camera,
				raw_scene.translations,
				raw_scene.scalings,
				raw_scene.rotations,
				recursion_depth,
				num_area_lights,
				this->light_sources.area_lights));
		}
	}

	planes = _planes;
	world = std::make_shared<BVH<TLASBox>>(objects);
}

Scene::~Scene() {}
