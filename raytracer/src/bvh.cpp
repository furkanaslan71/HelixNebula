#include "../include/bvh.h"
#include "../external/gsl/gsl"

struct TreeBVHNode {
  AABB bbox;
  TreeBVHNode* left = nullptr;
  TreeBVHNode* right = nullptr;
  int primitives_offset = 0;
  int primitive_count = 0;
  int split_axis = 0;
};

BVH::BVH(std::vector<std::shared_ptr<Hittable>>& _primitives)
  : primitives_(_primitives) 
{
  buildBVH();
}

TreeBVHNode* buildTreeRecursive(
  std::vector<std::shared_ptr<Hittable>>& objects,
  int start,
  int end)
{
  Expects(start < end);
  TreeBVHNode* node = new TreeBVHNode();

  int primitive_count = end - start;
  if (primitive_count == 1)
  {
    //leaf node
    node->bbox = objects[start]->getAABB();
    node->primitive_count = primitive_count;
    node->primitives_offset = start;
    return node;
  }

  AABB centroid_bbox;
  for (int i = start; i < end; i++)
    centroid_bbox.expand(objects[i]->getAABB().center());

  if (centroid_bbox.x.max == centroid_bbox.x.min &&
    centroid_bbox.y.max == centroid_bbox.y.min &&
    centroid_bbox.z.max == centroid_bbox.z.min)
  {
    // This means too many objects are located in one spot
    // put all of them together in a leaf node
    AABB bbox;
    for (int i = start; i < end; i++)
      bbox.expand(objects[i]->getAABB());

    node->bbox = bbox;
    node->primitives_offset = start;
    node->primitive_count = primitive_count;
    return node;
  }

  int longest_axis = centroid_bbox.longest_axis();
  int mid = (start + end) / 2;

  std::nth_element(&objects[start], &objects[mid], &objects[end - 1] + 1,
    [longest_axis](const std::shared_ptr<Hittable>& a,
      const std::shared_ptr<Hittable>& b) {
        return a->getAABB().center()[longest_axis] <
          b->getAABB().center()[longest_axis];
    });

  node->left = buildTreeRecursive(objects, start, mid);
  node->right = buildTreeRecursive(objects, mid, end);
  node->split_axis = longest_axis;
  node->primitive_count = 0;
  node->bbox = AABB(node->left->bbox, node->right->bbox);
  return node;
}

int BVH::buildFlatBVH(TreeBVHNode* node, int& offset)
{
  Expects(node);
  LinearBVHNode* linear_node = &linear_nodes_[offset];
  linear_node->bbox = node->bbox;
  int current_offset = offset++;

  if (node->primitive_count > 0)
  {
    linear_node->primitives_offset = node->primitives_offset;
    linear_node->primitive_count = node->primitive_count;
  }
  else
  {
    linear_node->primitive_count = 0;
    linear_node->axis = node->split_axis;
    buildFlatBVH(node->left, offset);
    linear_node->right_child_offset = buildFlatBVH(node->right, offset);
  }
  delete node;
  return current_offset;
}

void BVH::buildBVH()
{
  Expects(primitives_.size() > 0);
  size_t size = primitives_.size();
  TreeBVHNode* root = buildTreeRecursive(primitives_, 0, size);
  linear_nodes_.resize(size * 2);
  int offset = 0;
  buildFlatBVH(root, offset);
  linear_nodes_.resize(offset);
}

bool BVH::intersect(const Ray& ray, Interval ray_t, HitRecord& rec) const
{
  if (linear_nodes_.empty()) 
    return false;

  bool hit_anything = false;
  int current_node_index = 0;

  // Use a simpler stack pointer logic
  int stack_ptr = 0;
  int nodes_to_visit[INTERSECTION_STACK_SIZE];

  while (true)
  {
    const LinearBVHNode* node = &linear_nodes_[current_node_index];

    // IMPORTANT: Check against current closest hit (ray_t.max)
    // If we found a hit at t=5, we don't care about boxes starting at t=6
    if (node->bbox.hit(ray, ray_t))
    {
      if (node->primitive_count > 0) // Leaf
      {
        for (int i = 0; i < node->primitive_count; i++)
        {
          // Check intersection
          if (primitives_[node->primitives_offset + i]->hit(ray, ray_t, rec))
          {
            hit_anything = true;
            ray_t.max = rec.t; // CRITICAL: Shrink the ray so we ignore farther objects!
          }
        }
        if (stack_ptr == 0) break;
        current_node_index = nodes_to_visit[--stack_ptr];
      }
      else // Inner
      {
        // 1. Determine traversal order (Front-to-Back)
        bool dir_is_neg = ray.direction[node->axis] < 0;

        // 2. We want to visit the NEAR child immediately.
        //    So we push the FAR child onto the stack.
        int near_index, far_index;

        // Note: left child is always current + 1
        if (dir_is_neg)
        {
          near_index = node->right_child_offset;
          far_index = current_node_index + 1;
        }
        else
        {
          near_index = current_node_index + 1;
          far_index = node->right_child_offset;
        }

        nodes_to_visit[stack_ptr++] = far_index; // Push Far
        current_node_index = near_index;         // Visit Near
      }
    }
    else // Missed Box
    {
      if (stack_ptr == 0) 
        break;
      current_node_index = nodes_to_visit[--stack_ptr];
    }
  }
  return hit_anything;
}

AABB BVH::getAABB() const 
{ 
  return linear_nodes_[0].bbox; 
}

