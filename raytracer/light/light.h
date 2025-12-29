#ifndef LIGHT_H
#define LIGHT_H
#include <vector>

#include "math_core/math_core.h"
#include "core/color.h"

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

struct LightSources {
	std::vector<PointLight> point_lights;
	std::vector<AreaLight> area_lights;
	Color ambient_light;
};

#endif // LIGHT_H
