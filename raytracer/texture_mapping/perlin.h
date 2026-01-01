#ifndef PERLIN_H
#define PERLIN_H

#include <cmath>
#include <algorithm>  
#include <random>     
#include "glm_config.h"
#include "type_structs.h"


inline glm::vec3 gradients[16] = {
    glm::vec3(1, 1, 0),
    glm::vec3(-1, 1, 0),
    glm::vec3(1, -1, 0),
    glm::vec3(-1, -1, 0),
    glm::vec3(1, 0, 1),
    glm::vec3(-1, 0, 1),
    glm::vec3(1, 0, -1),
    glm::vec3(-1, 0, -1),
    glm::vec3(0, 1, 1),
    glm::vec3(0, -1, 1),
    glm::vec3(0, 1, -1),
    glm::vec3(0, -1, -1),
    glm::vec3(1, 1, 0),
    glm::vec3(-1, 1, 0),
    glm::vec3(0, -1, 1),
    glm::vec3(0, -1, -1)
};

inline int table[16] = {
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15
};



// Shuffle the array
static void initPerlin()
{
  static thread_local std::mt19937 g(std::random_device{}());
  std::shuffle(std::begin(table), std::end(table), g);
}


inline glm::vec3 getGradient(int i, int j, int k)
{
  int idx;
  idx = table[abs(k) % 16];
  idx = table[abs(j + idx) % 16];
  idx = table[abs(i + idx) % 16];
  return gradients[idx];
}

inline float interpolation_function(float x)
{
	float abs_x = std::abs(x);
	if (abs_x >= 1)
		return 0.0f;

	return -6.0f * abs_x * abs_x * abs_x * abs_x * abs_x
		+ 15.0f * abs_x * abs_x * abs_x * abs_x
		- 10 * abs_x * abs_x * abs_x
		+ 1.0f;
}

inline float calculate_c(int i, int j, int k, const glm::vec3& p)
{
  float dx = p.x - i;
  float dy = p.y - j;
  float dz = p.z - k;

  glm::vec3 d = glm::vec3(dx, dy, dz);
  glm::vec3 g = getGradient(i, j, k);

  float c = interpolation_function(dx) * interpolation_function(dy) * interpolation_function(dz) * glm::dot(g, d);
  return c;
}

inline float perlinNoise(const glm::vec3& point,
                         float noise_scale,
                         NoiseConversion noise_conversion,
                         int num_octaves)
{
  float s = 0.0f;
  const auto p = point * noise_scale;
  for (int l = 0; l < num_octaves; l++)
  {
    float pow2_l = std::pow(2, l);
    float inv_pow2_l = 1.0f / pow2_l;
    const auto scaled_p = p * pow2_l;

    int i, j, k;
    i = floor(scaled_p.x);
    j = floor(scaled_p.y);
    k = floor(scaled_p.z);

    float c = 0.0f;
    for (int di = 0; di <= 1; ++di)
    {
      for (int dj = 0; dj <= 1; ++dj)
      {
        for (int dk = 0; dk <= 1; ++dk)
        {
          int x = i + di;
          int y = j + dj;
          int z = k + dk;
          c += calculate_c(x, y, z, scaled_p);
        }
      }
    }
    
    s += c * inv_pow2_l;

  }

  if (noise_conversion == NoiseConversion::linear)
  {
    s = (s + 1) / 2.0f;
  }
  else
  {
    s = std::abs(s);
  }

  return s;
}

#endif // ! PERLIN_H
