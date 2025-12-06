#include "tlas_box.h"


TLASBox::TLASBox(
	int _local_index,
	int _world_index,
	const std::vector<std::shared_ptr<Hittable>>& _local_space_objects,
	const std::vector<std::optional<glm::mat4>>& _transform_matrices,
	const std::vector<int>& _material_ids,
	const std::vector<Vec3>& _motion_blur
)
	: 
	local_index(_local_index), 
	world_index(_world_index),
	local_space_objects(_local_space_objects),
	transform_matrices(_transform_matrices),
	material_ids(_material_ids),
	motion_blur(_motion_blur)
{
	const Hittable& base_object = *local_space_objects[local_index];
	AABB box = base_object.getAABB();
	const std::optional<glm::mat4>&
		composite_transformation_matrix = transform_matrices[world_index];
	if (composite_transformation_matrix.has_value())
		box = box.transformBox(composite_transformation_matrix.value());

	
	if (!motion_blur[world_index].isZero())
	{
		Vec3 motion_blur = this->motion_blur[world_index];
		glm::mat4 identity = glm::mat4(1.0);
		glm::mat4 translate = glm::translate(identity, motion_blur.toGlm3());
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

	const std::optional<glm::mat4>&
		composite_transformation_matrix = transform_matrices[world_index];

	const Hittable& base_object = *local_space_objects[local_index];

	Vec3 motion_offset(0, 0, 0);

	if (!motion_blur.empty() && !motion_blur[world_index].isZero())
	{
		motion_offset = motion_blur[world_index] * ray.time;
	}
	//motion_offset = Vec3(0.0, 0.0, 0.0);
	if (!composite_transformation_matrix.has_value())
	{
		Ray new_ray = ray;
		new_ray.origin = new_ray.origin - motion_offset;
		bool hit = base_object.hit(new_ray, ray_t, rec);
		if (hit)
		{
			rec.material_id = material_ids[world_index];
			rec.point = rec.point + motion_offset;
		}
			
		return hit;
	}

	const Vec3 world_origin = ray.origin - motion_offset;
	const Vec3 world_direction = ray.direction;


	const auto& matrix_val = composite_transformation_matrix.value();

	glm::mat4 inv_transform 
		= glm::inverse(matrix_val);

	glm::vec4 local_origin = inv_transform 
		* glm::vec4(world_origin.x, world_origin.y, world_origin.z, 1.0f);


	glm::vec4 local_direction = inv_transform 
		* glm::vec4(world_direction.x, world_direction.y, world_direction.z, 0.0f);

	Ray local_ray(
		Vec3(local_origin),
		Vec3(local_direction),
		ray.time
	);

	Interval local_ray_t(1e-8, INFINITY);


	if (!base_object.hit(local_ray, local_ray_t, rec))
		return false;

	glm::vec4 local_hit_point = glm::vec4(rec.point.x, rec.point.y, rec.point.z, 1.0f);
	glm::vec4 world_hit_point = matrix_val * local_hit_point;
	rec.point = Vec3(world_hit_point) + motion_offset;
	// Transform normal
	glm::vec4 local_normal = glm::vec4(rec.normal.x, rec.normal.y, rec.normal.z, 0.0f);
	glm::mat4 inv_transpose = glm::transpose(inv_transform);
	glm::vec4 world_normal = inv_transpose * local_normal;
	rec.normal = Vec3(world_normal).normalize();

	double local_t = rec.t;
	double world_t = (rec.point - ray.origin).dot(ray.direction) / (ray.direction.length() * ray.direction.length());


	if (world_t < ray_t.min || world_t > ray_t.max)
		return false;


	rec.t = world_t;

	rec.set_front_face(ray);
	rec.material_id = material_ids[world_index];
	return true;
}

AABB TLASBox::getAABB() const
{
	return bounding_box;
}