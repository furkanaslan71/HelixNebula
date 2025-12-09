#ifndef TLAS_BOX_H
#define TLAS_BOX_H
#include <optional>
#include <memory>
#include <vector>

#include "geometry.h"
#include "core/hittable.h"


class TLASBox{
public:
	TLASBox(int _world_index, 
					 std::vector<ObjectContext>* _object_contexes,
					Geometry* _geometry);

	bool hit(const Ray& ray, Interval ray_t, HitRecord& rec) const;
	AABB getAABB() const;

	 int world_index;
	
private:
	std::vector<ObjectContext>* object_contexes;
	AABB bounding_box;
	Geometry* geometry;
};

#endif
