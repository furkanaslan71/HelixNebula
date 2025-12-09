#include "ray.h"

Ray::Ray() {}

Ray::~Ray() {}

Ray::Ray(const  glm::vec3& _origin, const  glm::vec3& _direction, double _time)
  : origin(_origin), 
  direction(_direction),
	time(_time)
{
    //direction = glm::normalize(direction);
		inv_direction = 1.0f / direction;
		sign[0] = (inv_direction.x < 0);
		sign[1] = (inv_direction.y < 0);
		sign[2] = (inv_direction.z < 0);
}

void Ray::perturb(float roughness)
{
	glm::vec3 u, v;
	createONB(direction, u, v);
	float psi1 = generateRandomFloat(0, 1);
	float psi2 = generateRandomFloat(0, 1);
	//std::cout << psi1 << psi2 << std::endl;
	direction = direction + (u * (psi1 - 0.5f) + v * (psi2 - 0.5f)) * roughness;
	direction = glm::normalize(direction);
}
