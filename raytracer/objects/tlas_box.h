#ifndef TLAS_BOX_H
#define TLAS_BOX_H
#include "../include/hittable.h"
#include <optional>
#include <memory>

class TLASBox : public Hittable {
public:
	TLASBox(
		int _local_index,
		int _world_index, 
		const std::vector<std::shared_ptr<Hittable>>& _local_space_objects,
		const std::vector<std::optional<glm::mat4>>& _transform_matrices,
		const std::vector<int>& _material_ids,
		const std::vector<glm::vec3>& _motion_blur
	);
	const int local_index;
	const int world_index;
	bool hit(const Ray& ray, Interval ray_t, HitRecord& rec) const override;
	AABB getAABB() const override;
private:
	AABB bounding_box;
	const std::vector<std::shared_ptr<Hittable>>& local_space_objects;
	const std::vector<std::optional<glm::mat4>>& transform_matrices;
	const std::vector<int>& material_ids;
	const std::vector<glm::vec3>& motion_blur;
};

#endif
