#include "../include/aabb.h"
#include "../include/vec3.h"

AABB::AABB() {};
AABB::AABB(const glm::vec3& p1, const glm::vec3& p2)
{
    x = Interval(p1.x, p2.x);
    y = Interval(p1.y, p2.y);
    z = Interval(p1.z, p2.z);
    thicken();
}

AABB::AABB(const AABB& box1, const AABB& box2)
{
    x = Interval(box1.x, box2.x);
    y = Interval(box1.y, box2.y);
    z = Interval(box1.z, box2.z);
    thicken();
}


const void AABB::thicken()
{
    x.thicken();
    y.thicken();
    z.thicken();
}

const Interval& AABB::axis(int i) const
{
    if (i == 0) return x;
    if (i == 1) return y;
    else return z;
}

bool AABB::hit(const Ray& ray, Interval ray_t) const
{
    for (int i = 0; i < 3; i++)
    {
        double t0 = (axis(i).min - ray.origin[i]) / ray.direction[i];
        double t1 = (axis(i).max - ray.origin[i]) / ray.direction[i];
        if (t0 > t1) std::swap(t0, t1);

        // check whether overlaps
        if (ray_t.min < t0) ray_t.min = t0;
        if (ray_t.max > t1) ray_t.max = t1;

        if (ray_t.max <= ray_t.min) return false; // means no overlap between three axices.
    }
    return true;
}

const Interval& AABB::operator[](int axis) const
{
    switch (axis) {
        case 0:
            return x;
        case 1:
            return y;
        case 2:
            return z;
        default:
            throw std::out_of_range("Invalid axis index");
    }
}

AABB AABB::transformBox(const glm::mat4& matrix) const
{
  glm::vec3 corners[8];
    corners[0] = glm::vec3(x.min, y.min, z.min);
    corners[1] = glm::vec3(x.min, y.min, z.max);
    corners[2] = glm::vec3(x.min, y.max, z.min);
    corners[3] = glm::vec3(x.min, y.max, z.max);
    corners[4] = glm::vec3(x.max, y.min, z.min);
    corners[5] = glm::vec3(x.max, y.min, z.max);
    corners[6] = glm::vec3(x.max, y.max, z.min);
    corners[7] = glm::vec3(x.max, y.max, z.max);
    AABB transformedBox;
    for (auto& corner : corners)
    {
      corner = glm::vec3(matrix * glm::vec4(corner, 1.0));
    }

// The followng section is pure autism by erenjanje35
// let it stay for the fun of it

#define COMP(axis) [](const glm::vec3& a, const glm::vec3& b) {return a.axis < b.axis;}
#define FIND(kind, axis) double kind##_##axis = std::kind##_element(corners, corners + 8, COMP(axis))->axis
    FIND(min, x);
    FIND(max, x);
		FIND(min, y);
		FIND(max, y);
		FIND(min, z);
		FIND(max, z);
#undef COMP
#undef FIND
    transformedBox.x = Interval(min_x, max_x);
    transformedBox.y = Interval(min_y, max_y);
		transformedBox.z = Interval(min_z, max_z);

		return transformedBox;
}

void AABB::expand(const glm::vec3& p)
{
  x.min = std::min(x.min, (double)p.x);
  x.max = std::max(x.max, (double)p.x);
  y.min = std::min(y.min, (double)p.y);
  y.max = std::max(y.max, (double)p.y);
  z.min = std::min(z.min, (double)p.z);
  z.max = std::max(z.max, (double)p.z);
}

void AABB::expand(const AABB& other)
{
  x.min = std::min(x.min, other.x.min);
  x.max = std::max(x.max, other.x.max);
  y.min = std::min(y.min, other.y.min);
  y.max = std::max(y.max, other.y.max);
  z.min = std::min(z.min, other.z.min);
  z.max = std::max(z.max, other.z.max);
}

glm::vec3 AABB::center() const
{
  return glm::vec3(x.mid(), y.mid(), z.mid());
}

int AABB::longest_axis() const
{
  double dx = x.max - x.min;
  double dy = y.max - y.min;
  double dz = z.max - z.min;

  if (dx > dy && dx > dz) return 0;
  if (dy > dz) return 1;
  return 2;
}
