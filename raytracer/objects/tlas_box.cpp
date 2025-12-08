#include "tlas_box.h"
#define IS_ZERO (motion_blur[world_index].x == 0 && motion_blur[world_index].y == 0 && motion_blur[world_index].z == 0)

TLASBox::TLASBox(
	int _world_index,
	const std::vector<ObjectContext>& _object_contexes,
	std::shared_ptr<Hittable> _geometry
)
	: 
	world_index(_world_index),
	object_contexes(_object_contexes),
	geometry(_geometry)
{
	AABB box = geometry->getAABB();
	const auto& o_context = object_contexes[world_index];
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
	if (!bounding_box.hit(ray, ray_t))
		return false;

	const auto& o_context = object_contexes[world_index];
	const auto& composite_transformation_matrix = o_context.transform_matrix;

	const Hittable& base_object = *geometry;

	const auto& material_id = o_context.material_id;

	const auto& motion_blur = o_context.motion_blur;

	glm::vec3 motion_offset(0, 0, 0);

	if(!isZero(motion_blur))
		motion_offset = motion_blur * (float)ray.time;

	if (!composite_transformation_matrix.has_value())
	{
		Ray new_ray = ray;
		new_ray.origin = new_ray.origin - motion_offset;
		bool hit = base_object.hit(new_ray, ray_t, rec);
		if (hit)
		{
			rec.material_id = material_id;
			rec.point = rec.point + motion_offset;
		}
			
		return hit;
	}

	const glm::vec3 world_origin = ray.origin - motion_offset;
	const glm::vec3 world_direction = ray.direction;


	const auto& matrix_val = composite_transformation_matrix.value();

	glm::mat4 inv_transform = o_context.inverse_transform.value();

	glm::vec4 local_origin = inv_transform 
		* glm::vec4(world_origin.x, world_origin.y, world_origin.z, 1.0f);


	glm::vec4 local_direction = inv_transform 
		* glm::vec4(world_direction.x, world_direction.y, world_direction.z, 0.0f);

	Ray local_ray(
		glm::vec3(local_origin),
		glm::vec3(local_direction),
		ray.time
	);

	Interval local_ray_t(1e-8, INFINITY);


	if (!base_object.hit(local_ray, local_ray_t, rec))
		return false;

	glm::vec4 local_hit_point = glm::vec4(rec.point.x, rec.point.y, rec.point.z, 1.0f);
	glm::vec4 world_hit_point = matrix_val * local_hit_point;
	rec.point = glm::vec3(world_hit_point) + motion_offset;
	// Transform normal
	glm::vec4 local_normal = glm::vec4(rec.normal.x, rec.normal.y, rec.normal.z, 0.0f);
	glm::mat4 inv_transpose = glm::transpose(inv_transform);
	glm::vec4 world_normal = inv_transpose * local_normal;
	rec.normal = glm::normalize(glm::vec3(world_normal));

	double local_t = rec.t;
	double world_t = glm::dot((rec.point - ray.origin), ray.direction) / (glm::length(ray.direction) * glm::length(ray.direction));


	if (world_t < ray_t.min || world_t > ray_t.max)
		return false;


	rec.t = world_t;

	rec.set_front_face(ray);
	rec.material_id = material_id;
	return true;
}

AABB TLASBox::getAABB() const
{
	return bounding_box;
}