#ifndef HITTABLE_H
#define HITTABLE_H
#include <optional>

#include "glm_config.h"
#include "core/ray.h"

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

struct ObjectContext {
  std::optional<glm::mat4> transform_matrix;
  std::optional<glm::mat4> inverse_transform;
  glm::vec3 motion_blur;
  int material_id;
};


//class Hittable {
//public:
//  virtual ~Hittable() = default;
//  virtual bool hit(const Ray& ray, Interval ray_t, HitRecord& rec) const = 0;
//  virtual AABB getAABB() const = 0;
//};

#endif //HITTABLE_H
