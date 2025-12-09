#ifndef PLANE_H
#define PLANE_H
#include <string>
#include <optional>

#include "parser/parser.hpp"
#include "math_core/math_core.h"
#include "core/ray.h"
#include "core/hittable.h"
#include "core/interval.h"


class Plane {
public:
	Plane();
	Plane(
		const Plane_& _plane,
		const std::vector<glm::vec3>& _vertex_data,
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