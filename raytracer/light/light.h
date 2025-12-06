#ifndef LIGHT_H
#define LIGHT_H

#include <vector>
#include "../include/vec3.h"
#include "../include/color.h"
#include "../random/my_random.h"


struct PointLight {
	int id;
	Vec3 position;
	Color intensity;
};

struct AreaLight {
	int id;
	Vec3 position;
	Vec3 normal;
	float edge;
	Vec3 radiance;
	Vec3 u, v;
	void generateApertureSamples(
		std::vector<Vec3>& out_samples, int num_samples
	) const
	{
		std::vector<std::pair<float, float>> samples
			= generateJitteredSamples(num_samples);

		out_samples.clear();
		for (const auto& [x, y] : samples)
		{
			Vec3 area_light_sample
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
