#include "plane.h"

Plane::Plane()
		: id(-1), material_id(-1), point(glm::vec3(0, 0, 0)), normal(glm::vec3(0, 1, 0))
	{
}

Plane::Plane(
	const Plane_& _plane,
	const std::vector<Vec3f_>& _vertex_data,
	glm::vec3 _motion_blur
)
	: id(_plane.id), material_id(_plane.material_id),
	point(glm::vec3(_vertex_data[_plane.point_vertex_id])),
	normal(glm::normalize(glm::vec3(_plane.normal))),
	motion_blur(_motion_blur)
{
	composite_transformation_matrix = _plane.transform_matrix;
	if (composite_transformation_matrix.has_value())
	{
		//transform point and normal
		glm::vec4 glm_point = glm::vec4(static_cast<float>(point.x),
			static_cast<float>(point.y),
			static_cast<float>(point.z),
			1.0f);
		glm::vec4 transformed_point = composite_transformation_matrix.value() * glm_point;
		point = glm::vec3(transformed_point);
		glm::vec4 glm_normal = glm::vec4(static_cast<float>(normal.x),
			static_cast<float>(normal.y),
			static_cast<float>(normal.z),
			0.0f);
		glm::vec4 transformed_normal = glm::transpose(glm::inverse(composite_transformation_matrix.value())) * glm_normal;
		normal = glm::normalize(glm::vec3(transformed_normal));
	}
	
}
bool Plane::hit(const Ray& ray, const Interval& interval, HitRecord& rec) const
{
	double denom = glm::dot(normal, ray.direction);
	if (std::abs(denom) > 1e-6)
	{
		glm::vec3 p0l0 = point - ray.origin;
		double t = glm::dot(p0l0, normal) / denom;
		if (t >= interval.min && t <= interval.max)
		{
			rec.t = t;
			rec.point = ray.origin + ray.direction * (float)t;
			rec.normal = normal;
			rec.material_id = material_id;
			rec.set_front_face(ray);
			return true;
		}
	}
	return false;
}