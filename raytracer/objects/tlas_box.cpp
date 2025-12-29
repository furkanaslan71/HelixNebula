#include "tlas_box.h"
#define IS_ZERO (motion_blur[world_index].x == 0 && motion_blur[world_index].y == 0 && motion_blur[world_index].z == 0)

TLASBox::TLASBox(
	int _world_index,
	 std::vector<ObjectContext>* _object_contexes,
	Geometry* _geometry
)
	: 
	world_index(_world_index),
	object_contexes(_object_contexes),
	geometry(_geometry)
{
	AABB box = geometry->getAABB();
	const auto& o_context = (*object_contexes)[world_index];
	const std::optional<glm::mat4>& composite_transformation_matrix = o_context.transform_matrix;
	if (composite_transformation_matrix.has_value())
		box = box.transformBox(composite_transformation_matrix.value());


	const auto& motion_blur = o_context.motion_blur;
	if (isZero(motion_blur))
	{
		glm::mat4 identity = glm::mat4(1.0);
		glm::mat4 translate = glm::translate(identity, motion_blur);
		AABB t_0_pos = box;
		AABB t_1_pos = t_0_pos.transformBox(translate);
		box.expand(t_1_pos);
	}

	bounding_box = box;
}

bool TLASBox::hit(const Ray& ray, Interval ray_t, HitRecord& rec) const
{
	const auto& o_context = (*object_contexes)[world_index];
	const auto& composite_transformation_matrix = o_context.transform_matrix;

	glm::vec3 motion_offset = o_context.motion_blur * (float)ray.time;

	if (!composite_transformation_matrix.has_value())
	{
		Ray new_ray(ray.origin - motion_offset, ray.direction, ray.time);
		bool hit = geometry->hit(new_ray, ray_t, rec);
		if (hit)
		{
			rec.material_id = o_context.material_id;
			rec.point += motion_offset;
		}
			
		return hit;
	}

	const glm::mat4& inv_transform = o_context.inverse_transform.value();

	glm::vec4 local_origin = inv_transform * glm::vec4(ray.origin - motion_offset, 1.0f);
	glm::vec4 local_direction = inv_transform * glm::vec4(ray.direction, 0.0f);

	Ray local_ray(glm::vec3(local_origin), glm::vec3(local_direction),ray.time);


	if (!geometry->hit(local_ray, ray_t, rec))
		return false;

	const auto& matrix_val = composite_transformation_matrix.value();

	rec.point = glm::vec3(composite_transformation_matrix.value() * glm::vec4(rec.point, 1.0f));
	rec.point += motion_offset;

	glm::mat4 inv_transpose = glm::transpose(inv_transform);
	rec.normal = glm::normalize(glm::vec3(inv_transpose * glm::vec4(rec.normal, 0.0f)));

	rec.set_front_face(ray);
	rec.material_id = o_context.material_id;
	return true;
}

AABB TLASBox::getAABB() const
{
	return bounding_box;
}