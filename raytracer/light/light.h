#ifndef LIGHT_H
#define LIGHT_H
#include <vector>

#include "math_core/math_core.h"
#include "core/color.h"

enum class EnvMapType {
	latlong,
	probe
};

enum class Sampler {
	uniform,
	cosine
};

struct PointLight {
	int id;
	glm::vec3 position;
	Color intensity;
};

struct AreaLight {
	int id;
	glm::vec3 position;
	glm::vec3 normal;
	float edge;
	glm::vec3 radiance;
	glm::vec3 u, v;

	void generateAreaLightSamples(
		std::vector<glm::vec3>& out_samples, int num_samples
	) const
	{
		std::vector<std::pair<float, float>> samples
			= generateJitteredSamples(num_samples);

		out_samples.clear();
		for (const auto& [x, y] : samples)
		{
			glm::vec3 area_light_sample
				= position + (u * (x - 0.5f) + v * (y - 0.5f)) * edge;
			out_samples.emplace_back(area_light_sample);
		}
	}
};

struct DirectionalLight {
	DirectionalLight() = default;
	explicit DirectionalLight(const DirectionalLight_& _dl)
		: direction(_dl.direction), radiance(_dl.radiance), id(_dl.id){}

	glm::vec3 direction;
	glm::vec3 radiance;
	int id;
};

struct SpotLight {
	SpotLight()= default;
	SpotLight(const SpotLight_& _sl)
		:
	position(_sl.position),
	direction(_sl.direction),
	intensity(_sl.intensity),
	coverage_angle(_sl.coverage_angle),
	falloff_angle(_sl.falloff_angle),
	id(_sl.id)
	{};

	glm::vec3 position;
	glm::vec3 direction;
	glm::vec3 intensity;
	float coverage_angle;
	float falloff_angle;
	int id;
};

struct EnvironmentLight {
	EnvironmentLight() = default;
	EnvironmentLight(const SphericalDirectionalLight_& _sdl, Image* _img)
		: id(_sdl.id), img(_img)
	{
		if (_sdl.sampler == "uniform")
			sampler = Sampler::uniform;
		else if (_sdl.sampler == "cosine")
			sampler = Sampler::cosine;
		else
			throw std::runtime_error("Unknown sampler type");

		if (_sdl.type == "latlong")
			type = EnvMapType::latlong;
		else if (_sdl.type == "probe")
			type = EnvMapType::probe;
		else
			throw std::runtime_error("Unknown env_map type type");
	}
 	Image* img = nullptr;
	EnvMapType type;
	Sampler sampler;
	int id;
};

struct LightSources {
	std::vector<PointLight> point_lights;
	std::vector<AreaLight> area_lights;
	std::vector<DirectionalLight> directional_lights;
	std::vector<SpotLight> spot_lights;
	EnvironmentLight env_light;
	Color ambient_light;
};

#endif // LIGHT_H
