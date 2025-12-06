#ifndef RAY_H
#define RAY_H

#include "vec3.h"
#include "../random/my_random.h"

class Ray {
public:
    Vec3 origin;
    Vec3 direction;
    double time;
    Ray();
    ~Ray();
    Ray(const Vec3& _origin, const Vec3& _direction, double _time);

    void perturb(float roughness);
};
#endif //RAY_H
