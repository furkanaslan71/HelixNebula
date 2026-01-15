#ifndef HITTABLE_H
#define HITTABLE_H

#include <vector>

#include "glm_config.h"
#include "core/ray.h"
#include "material/material.h"
#include "texture_mapping/texture_data.h"

struct TangentVectors {
  glm::vec3 u;
  glm::vec3 v;
};

typedef struct HitRecord{
  TangentVectors surface_tangents;
  glm::vec3 point;
  glm::vec3 normal;
  Material* material;
  std::vector<Texture*> textures;
  glm::vec2 uv;
  double t;
	bool front_face;
  std::optional<glm::vec3> radiance;
  
  void set_front_face(const Ray& r)
  {
    front_face = glm::dot(r.direction, normal) < 0;
  }
}HitRecord;



#endif //HITTABLE_H
