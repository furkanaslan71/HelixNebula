#ifndef BVH_H
#define BVH_H
#include <memory>
#include <vector>
#include <algorithm>

#include "core/hittable.h"
#include "objects/geometry_concepts.h"
#include "external/gsl/gsl"

struct alignas(32) LinearBVHNode {
	AABB bbox;
	union {
		int primitives_offset;
		int right_child_offset;
	};
	uint16_t primitive_count;
	uint8_t axis;
	uint8_t padding;
};

struct TreeBVHNode;

template <GeometryConcept T>
class BVH{
public:
	BVH(std::vector<T>& _primitives);

	bool intersect(const Ray& ray, Interval ray_t, HitRecord& rec) const;
	AABB getAABB() const;
	void buildBVH();
private:
	
	int buildFlatBVH(TreeBVHNode* node, int& offset);

	std::vector<T>& primitives_;
	std::vector<LinearBVHNode> linear_nodes_;
};

#include "bvh.tpp"

#endif // !BVH_H