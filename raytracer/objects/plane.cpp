#include "plane.h"

Plane::Plane()
		: id(-1), material_id(-1), point(Vec3(0, 0, 0)), normal(Vec3(0, 1, 0))
	{
}

Plane::Plane(
	const Plane_& _plane,
	const std::vector<Vec3f_>& _vertex_data,
	Vec3 _motion_blur
)
	: id(_plane.id), material_id(_plane.material_id),
	point(Vec3(_vertex_data[_plane.point_vertex_id])),
	normal(Vec3(_plane.normal).normalize()),
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
		point = Vec3(
			static_cast<double>(transformed_point.x),
			static_cast<double>(transformed_point.y),
			static_cast<double>(transformed_point.z)
		);
		glm::vec4 glm_normal = glm::vec4(static_cast<float>(normal.x),
			static_cast<float>(normal.y),
			static_cast<float>(normal.z),
			0.0f);
		glm::vec4 transformed_normal = glm::transpose(glm::inverse(composite_transformation_matrix.value())) * glm_normal;
		normal = Vec3(
			static_cast<double>(transformed_normal.x),
			static_cast<double>(transformed_normal.y),
			static_cast<double>(transformed_normal.z)
		).normalize();
	}
	
}
bool Plane::hit(const Ray& ray, const Interval& interval, HitRecord& rec) const
{
	double denom = normal.dot(ray.direction);
	if (std::abs(denom) > 1e-6)
	{
		Vec3 p0l0 = point - ray.origin;
		double t = p0l0.dot(normal) / denom;
		if (t >= interval.min && t <= interval.max)
		{
			rec.t = t;
			rec.point = ray.origin + ray.direction * t;
			rec.normal = normal;
			rec.material_id = material_id;
			rec.set_front_face(ray);
			return true;
		}
	}
	return false;
}