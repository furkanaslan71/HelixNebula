#ifndef BVH_H
#define BVH_H

#include "hittable.h"
#include <memory>
#include <vector>
#include <algorithm>
#include <random>

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

class BVH{
public:
	BVH(std::vector<std::shared_ptr<Hittable>>& _primitives);

	bool intersect(const Ray& ray, Interval ray_t, HitRecord& rec) const;
	AABB getAABB() const;
	void buildBVH();
private:
	
	int buildFlatBVH(TreeBVHNode* node, int& offset);

	std::vector<std::shared_ptr<Hittable>>& primitives_;
	std::vector<LinearBVHNode> linear_nodes_;
};

#endif // !BVH_H