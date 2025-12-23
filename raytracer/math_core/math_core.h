#ifndef MATH_CORE_H
#define MATH_CORE_H

constexpr auto INV_255f = 1.0f/255.0f;

#include "glm_config.h"
#include "vector_utility.h"
#include "transformations/transformations.h"
#include "random/my_random.h"

struct BackgroundCameraData {
  glm::vec3 cam_forward;
  glm::vec3 cam_right;
  glm::vec3 cam_up;

  float tan_half_fov_x;
  float tan_half_fov_y;

};

#endif //!MATH_CORE_H