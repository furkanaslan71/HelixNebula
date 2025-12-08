#ifndef TLAS_BOX_H
#define TLAS_BOX_H
#include "../include/hittable.h"
#include <optional>
#include <memory>

class TLASBox : public Hittable {
public:
	TLASBox(int _world_index, 
					const std::vector<ObjectContext>& _object_contexes,
					std::shared_ptr<Hittable> _geometry);

	bool hit(const Ray& ray, Interval ray_t, HitRecord& rec) const override;
	AABB getAABB() const override;

	const int world_index;
	
private:
	const std::vector<ObjectContext>& object_contexes;
	AABB bounding_box;
	const std::shared_ptr<Hittable> geometry;
};

#endif
