#ifndef MESH_H
#define MESH_H

#include "core/hittable.h"
#include "parser/parser.hpp"
#include "acceleration/bvh.h"
#include "triangle.h"


template<Shading mode, TextureLookup tex>
class Mesh {
public:
	Mesh(int _id, const std::vector<Triangle_>& _faces,
			 const std::vector<glm::vec3>& vertex_data, const std::vector<glm::vec2> uv_data);

	bool hit(const Ray& ray, Interval ray_t, HitRecord& rec) const;

	AABB getAABB() const;

	int id;
private:
	std::vector<TriangleNew<mode, tex>> faces;
	AABB bounding_box;
	std::shared_ptr<BVH<TriangleNew<mode, tex>>> bvh;
};

#include "mesh.tpp"


#endif