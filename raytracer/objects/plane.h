#ifndef PLANE_H
#define PLANE_H


#include "../include/parser.hpp"
#include "../math_core/math_core.h"
#include "../include/ray.h"
#include "../include/hittable.h"
#include "../include/interval.h"
#include <optional>

class Plane {
public:
	Plane();
	Plane(
		const Plane_& _plane,
		const std::vector<Vec3f_>& _vertex_data,
		glm::vec3 _motion_blur
	);
	bool hit(const Ray& ray, const Interval& interval, HitRecord& rec) const;
	int id;
	int material_id;
	glm::vec3 point;
	glm::vec3 normal;
	glm::vec3 motion_blur;
	std::optional<glm::mat4> composite_transformation_matrix;
};
#endif // PLANE_H