#include "distribution_camera.h"

void DistributionCamera::generateApertureSamples(
	std::vector<glm::vec3>& out_samples) const
{
	std::vector<std::pair<float, float>> samples
		= generateJitteredSamples(num_samples);

	out_samples.clear();
	for (const auto& [x, y] : samples)
	{
		glm::vec3 aperture_sample = position + (u * (x - 0.5f) + v * (y - 0.5f)) * aperture_size;
		out_samples.emplace_back(aperture_sample);
	}
}

glm::vec3 DistributionCamera::calculateDir(
	const glm::vec3& pixel_sample, const glm::vec3& a) const
{
	glm::vec3 e = position;
	glm::vec3 dir = pixel_sample - e;
	float t_fp = focus_distance / glm::dot(dir, -w);
	glm::vec3 p = e + dir * t_fp;
	glm::vec3 d = p - a;
	return d;
}

void DistributionCamera::render(
	const BaseRayTracer& rendering_technique,
	std::vector<std::vector<Color>>& image) const

{
	image.resize(image_height, std::vector<Color>(image_width, Color(0, 0, 0)));

#if DEBUG_PIXEL
	Vec3 pixel_center = q + su * (WIDTH + 0.5) + sv * (HEIGHT + 0.5);
	Ray primary_ray(position, (pixel_center - position).normalize());

	Color pixel_color = rendering_technique.traceRay(primary_ray);
	pixel_color = pixel_color.clamp();
	image[HEIGHT][WIDTH] = pixel_color;
#else

#if !MULTI_THREADING
	for (int i = 0; i < image_height; ++i)
	{
		for (int j = 0; j < image_width; ++j)
		{
			Vec3 pixel_center = q + su * (j + 0.5) + sv * (i + 0.5);
			Ray primary_ray(position, (pixel_center - position).normalize());

			Color pixel_color = rendering_technique.traceRay(primary_ray);
			pixel_color = pixel_color.clamp();
			image[i][j] = pixel_color;
		}
	}
#else
	const int numThreads = std::thread::hardware_concurrency();
	std::vector<std::thread> threads(numThreads);
	for (int threadId = 0; threadId < numThreads; threadId++)
	{
		threads[threadId] = std::thread([this, threadId, numThreads,
			&rendering_technique, &image]() {
				std::vector<glm::vec3> pixel_samples;
				pixel_samples.reserve(num_samples);

				std::vector<glm::vec3> aperture_samples;
				aperture_samples.reserve(num_samples);

				std::vector<std::vector<std::vector<glm::vec3>>> area_light_samples;
				area_light_samples.resize(recursion_depth + 1);
				for (int l = 0; l < recursion_depth + 1; l++)
				{
					area_light_samples[l].resize(num_area_lights);
					for (int m = 0; m < num_area_lights; m++)
					{
						area_light_samples[l][m].reserve(num_samples);
					}
				}

				std::mt19937 rng(std::random_device{}());
				for (int i = threadId; i < image_height; i += numThreads)
				{
					for (int j = 0; j < image_width; ++j)
					{
						Color pixel_color = Color(0.0, 0.0, 0.0);
						generatePixelSamples(i, j, pixel_samples);
						generateApertureSamples(aperture_samples);
						std::shuffle(aperture_samples.begin(), aperture_samples.end(), rng);

						for (auto& depth : area_light_samples)
						{
							for (int f = 0; f < num_area_lights; f++)
							{
								area_lights[f].generateApertureSamples(depth[f], num_samples);
								std::shuffle(depth[f].begin(), depth[f].end(), rng);
							}
						}

						for (int k = 0; k < num_samples; k++)
						{
							glm::vec3 a = (aperture_size > 0.0) ? aperture_samples[k] : position;
							glm::vec3 dir = glm::normalize(calculateDir(pixel_samples[k], a));

							Ray primary_ray(a, dir, generateRandomFloat(0, 1));

							RenderContext context;
							context.area_light_samples = &area_light_samples; // Point to our data
							context.sample_index = k; // Tell the shader we are on sample 'k'

							pixel_color += rendering_technique.traceRay(primary_ray, context);
						}
						pixel_color = pixel_color / (float)num_samples;
						pixel_color = pixel_color.clamp();
						image[i][j] = pixel_color;
					}
				}
			});
	}
	for (int threadId = 0; threadId < numThreads; ++threadId)
	{
		threads[threadId].join();
	}
#endif

#endif

}