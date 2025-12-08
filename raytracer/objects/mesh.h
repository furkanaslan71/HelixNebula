#ifndef MESH_H
#define MESH_H
#include "../include/hittable.h"
#include "../include/parser.hpp"
#include "../include/bvh.h"
#include "triangle.h"


class Mesh {
public:
	Mesh(int _id, bool _smooth_shading, const std::vector<Triangle_>& _faces,
			 const std::vector<glm::vec3>& vertex_data);

	bool hit(const Ray& ray, Interval ray_t, HitRecord& rec) const;

	AABB getAABB() const;

	int id;
	bool smooth_shading;
private:
	std::vector<Triangle> faces;
	const AABB bounding_box;
	std::shared_ptr<BVH<Triangle>> bvh;
};

#endif