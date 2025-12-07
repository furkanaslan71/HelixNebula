#ifndef HITTABLE_H
#define HITTABLE_H

#include "../include/ray.h"
#include "interval.h"
#include "aabb.h"

typedef struct HitRecord{
  glm::vec3 point;
  glm::vec3 normal;
	bool front_face;
  int material_id;
  double t;
  void set_front_face(const Ray& r)
  {
    front_face = glm::dot(r.direction, normal) < 0;
  }
}HitRecord;

class Hittable {
public:
  virtual ~Hittable() = default;
  virtual bool hit(const Ray& ray, Interval ray_t, HitRecord& rec) const = 0;
  virtual AABB getAABB() const = 0;
};

#endif //HITTABLE_H
