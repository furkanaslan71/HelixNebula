#include "tlas_box.h"
#define IS_ZERO (motion_blur[world_index].x == 0 && motion_blur[world_index].y == 0 && motion_blur[world_index].z == 0)

template bool TLASBox::hit<true>(const Ray&, Interval, HitRecord&) const;
template bool TLASBox::hit<false>(const Ray&, Interval, HitRecord&) const;

TLASBox::TLASBox(Geometry* _geometry,
	Material* _material,
	std::optional<Transformation>* _transformation,
	glm::vec3 _motion_blur,
	std::vector<Texture*> _textures)
		:
	geometry(_geometry),
	material(_material),
	transformation(_transformation),
	motion_blur(_motion_blur)
{
	for (const auto tex_ptr : _textures)
	{
		Expects(tex_ptr);
		textures.emplace_back(tex_ptr);
	}
	AABB box = geometry->getAABB();

	if (transformation->has_value())
		box = box.transformBox(transformation->value().transform_matrix);

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

template <bool occlusion_only>
bool TLASBox::hit(const Ray& ray, Interval ray_t, HitRecord& rec) const
{
	glm::vec3 motion_offset = motion_blur * static_cast<float>(ray.time);

	if (!transformation->has_value())
	{
		Ray new_ray(ray.origin - motion_offset, ray.direction, ray.time);
		bool hit = geometry->hit<occlusion_only>(new_ray, ray_t, rec);
		if (hit)
		{
			if constexpr (occlusion_only) return true;
			rec.material = material;
			rec.point += motion_offset;
			rec.textures = this->textures;
		}
		return hit;
	}

	const glm::mat4& inv_transform = transformation->value().inverse_transform;

	glm::vec4 local_origin = inv_transform * glm::vec4(ray.origin - motion_offset, 1.0f);
	glm::vec4 local_direction = inv_transform * glm::vec4(ray.direction, 0.0f);

	Ray local_ray(glm::vec3(local_origin), glm::vec3(local_direction),ray.time);


	if (!geometry->hit<occlusion_only>(local_ray, ray_t, rec))
		return false;
	if constexpr (occlusion_only) return true;

	const auto& matrix_val = transformation->value().transform_matrix;

	rec.point = glm::vec3(matrix_val * glm::vec4(rec.point, 1.0f));
	rec.point += motion_offset;

	glm::mat4 inv_transpose = glm::transpose(inv_transform);
	rec.normal = glm::normalize(glm::vec3(inv_transpose * glm::vec4(rec.normal, 0.0f)));

	rec.set_front_face(ray);
	rec.material = material;
	rec.textures = this->textures;
	return true;
}

AABB TLASBox::getAABB() const
{
	return bounding_box;
}