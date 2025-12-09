#ifndef AABB_H
#define AABB_H

#include "interval.h"
#include "math_core/math_core.h"
#include "ray.h"

class AABB {
public:
    Interval x, y, z;

    AABB();
    AABB(const glm::vec3& p1, const glm::vec3& p2);
    AABB(const AABB& box1, const AABB& box2);

    const void thicken();
    const Interval& axis(int i) const;
    bool hit(const Ray& ray, Interval ray_t) const;
    const Interval& operator[](int axis) const;
		AABB transformBox(const glm::mat4& matrix) const;
    void expand(const glm::vec3& p);
    void expand(const AABB& other);
    glm::vec3 center() const;
    int longest_axis() const;
};

#endif // AABB_H