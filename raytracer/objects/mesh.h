#ifndef MESH_H
#define MESH_H

#include "core/hittable.h"
#include "parser/parser.hpp"
#include "acceleration/bvh.h"
#include "triangle.h"

class Mesh {
public:
	Mesh(int _id, const std::vector<Triangle_>& _faces,
			 const std::vector<glm::vec3>& vertex_data,
			 int vertex_offset, int texture_offset, bool is_smooth,
			 const std::vector<glm::vec2>& tex_data);

	template <bool occlusion_only>
	bool hit(const Ray& ray, Interval ray_t, HitRecord& rec) const;

	AABB getAABB() const;

	int id;
private:
	std::vector<Triangle> faces;
	AABB bounding_box;
	std::shared_ptr<BVH<Triangle>> bvh;
};


#endif
