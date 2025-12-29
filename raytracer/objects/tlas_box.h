#ifndef TLAS_BOX_H
#define TLAS_BOX_H
#include <optional>
#include <memory>
#include <vector>
#include <material/material.h>

#include "geometry.h"
#include "core/hittable.h"

struct Transformation {
	glm::mat4 transform_matrix;
	glm::mat4 inverse_transform;
};

class TLASBox{
public:
	TLASBox(Geometry* _geometry, Material* _material, std::optional<Transformation>* _transformation, glm::vec3 _motion_blur);
	bool hit(const Ray& ray, Interval ray_t, HitRecord& rec) const;
	AABB getAABB() const;

private:
	AABB bounding_box;
	glm::vec3 motion_blur;
	Geometry* geometry;
	Material* material;
	std::optional<Transformation>* transformation;
};

#endif
