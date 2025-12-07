#ifndef RAY_H
#define RAY_H

#include "vec3.h"
#include "../random/my_random.h"

class Ray {
public:
    glm::vec3 origin;
    glm::vec3 direction;
    double time;
    Ray();
    ~Ray();
    Ray(const  glm::vec3& _origin, const  glm::vec3& _direction, double _time);

    void perturb(float roughness);
};
#endif //RAY_H
