#include "scene.h"

#include "../include/ray.h"


Scene::Scene() {}

Scene::Scene(const Scene_& raw_scene, std::vector<std::shared_ptr<Hittable>>& objects)
	: background_color(raw_scene.background_color.x, raw_scene.background_color.y, raw_scene.background_color.z)
{
	for (const auto& raw_light : raw_scene.point_lights) 
	{
		auto transformation_matrix = calculateCompositeTransformationMatrix(
			raw_scene.translations,
			raw_scene.scalings,
			raw_scene.rotations,
			raw_light);
		glm::vec4 position_v4 = transformation_matrix * 
			glm::vec4(raw_light.position.x, 
				raw_light.position.y, 
				raw_light.position.z, 1.0f);
		Vec3 transformed_position = 
			Vec3(position_v4.x, position_v4.y, position_v4.z);

		light_sources.point_lights.push_back(PointLight(raw_light.id,
			Vec3(transformed_position),
			Color(raw_light.intensity.x, raw_light.intensity.y, raw_light.intensity.z)));
	}
	for (const auto& raw_light : raw_scene.area_lights)
	{
		auto transformation_matrix = calculateCompositeTransformationMatrix(
			raw_scene.translations,
			raw_scene.scalings,
			raw_scene.rotations,
			raw_light);
		glm::vec4 position_v4 = transformation_matrix *
			glm::vec4(raw_light.position.x,
				raw_light.position.y,
				raw_light.position.z, 1.0f);
		Vec3 transformed_position =
			Vec3(position_v4.x, position_v4.y, position_v4.z);

		Vec3 u, v;

		Vec3::createONB(raw_light.normal, u, v);

		light_sources.area_lights.push_back(
			AreaLight(
				raw_light.id,
				Vec3(transformed_position),
				Vec3(raw_light.normal),
				raw_light.edge,
				Vec3(raw_light.radiance.x, raw_light.radiance.y, raw_light.radiance.z),
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
			pinhole_cameras.emplace_back(raw_camera,
				raw_scene.translations,
				raw_scene.scalings,
				raw_scene.rotations,
				recursion_depth,
				num_area_lights,
				this->light_sources.area_lights
				);
		}
		else
		{
			distribution_cameras.emplace_back(raw_camera,
				raw_scene.translations,
				raw_scene.scalings,
				raw_scene.rotations,
				recursion_depth,
				num_area_lights,
				this->light_sources.area_lights);
		}
	}
		
	world = std::make_shared<BVH>(objects);
}

Scene::~Scene() {}

glm::mat4 Scene::calculateCompositeTransformationMatrix(
	const std::vector<Translation_>& translations,
	const std::vector<Scaling_>& scalings,
	const std::vector<Rotation_>& rotations,
	const PointLight_& raw_light)
{
	auto transformations = raw_light.transformations;
	auto res = glm::mat4(1.0f);
	for (const auto& transform : transformations)
	{
		auto identity = glm::mat4(1.0f);
		switch (transform[0])
		{
		case 't': // translation
		{
			int index = std::stoi(transform.substr(1));
			auto& t = translations[index];
			glm::vec3 translation = glm::vec3(t.tx, t.ty, t.tz);
			res = glm::translate(identity, translation) * res;
			break;
		}
		case 's': // scaling
		{
			int index = std::stoi(transform.substr(1));
			auto& s = scalings[index];
			glm::vec3 scaling = glm::vec3(s.sx, s.sy, s.sz);
			res = glm::scale(identity, scaling) * res;
			break;
		}
		case 'r': // rotation
		{
			int index = std::stoi(transform.substr(1));
			auto& r = rotations[index];
			glm::vec3 axis = glm::vec3(r.axis_x, r.axis_y, r.axis_z);
			res = glm::rotate(identity, glm::radians(r.angle), axis) * res;
			break;
		}
		default:
			break;
		}
	}
	return res;
}

glm::mat4 Scene::calculateCompositeTransformationMatrix(
	const std::vector<Translation_>& translations,
	const std::vector<Scaling_>& scalings,
	const std::vector<Rotation_>& rotations,
	const AreaLight_& raw_light)
{
	auto transformations = raw_light.transformations;
	auto res = glm::mat4(1.0f);
	for (const auto& transform : transformations)
	{
		auto identity = glm::mat4(1.0f);
		switch (transform[0])
		{
		case 't': // translation
		{
			int index = std::stoi(transform.substr(1));
			auto& t = translations[index];
			glm::vec3 translation = glm::vec3(t.tx, t.ty, t.tz);
			res = glm::translate(identity, translation) * res;
			break;
		}
		case 's': // scaling
		{
			int index = std::stoi(transform.substr(1));
			auto& s = scalings[index];
			glm::vec3 scaling = glm::vec3(s.sx, s.sy, s.sz);
			res = glm::scale(identity, scaling) * res;
			break;
		}
		case 'r': // rotation
		{
			int index = std::stoi(transform.substr(1));
			auto& r = rotations[index];
			glm::vec3 axis = glm::vec3(r.axis_x, r.axis_y, r.axis_z);
			res = glm::rotate(identity, glm::radians(r.angle), axis) * res;
			break;
		}
		default:
			break;
		}
	}
	return res;
}

