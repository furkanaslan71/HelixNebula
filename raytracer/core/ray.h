#ifndef RAY_H
#define RAY_H

#include "math_core/math_core.h"

class Ray {
public:
    glm::vec3 origin;
    glm::vec3 direction;
    glm::vec3 inv_direction;
    double time;
    int8_t sign[3];
    int8_t inside = 0;
    Ray();
    ~Ray();
    Ray(const  glm::vec3& _origin, const  glm::vec3& _direction, double _time);

    void perturb(float roughness);
};
#endif //RAY_H
