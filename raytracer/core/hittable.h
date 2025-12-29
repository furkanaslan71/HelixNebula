#ifndef HITTABLE_H
#define HITTABLE_H
#include <optional>
#include <vector>

#include "glm_config.h"
#include "core/ray.h"
#include "material/material.h"

typedef struct HitRecord{
  glm::vec3 point;
  glm::vec3 normal;
  double t;
  Material* material;
	bool front_face;
  
  void set_front_face(const Ray& r)
  {
    front_face = glm::dot(r.direction, normal) < 0;
  }
}HitRecord;



#endif //HITTABLE_H
