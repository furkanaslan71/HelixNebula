#include "../include/ray.h"

Ray::Ray() {}

Ray::~Ray() {}

Ray::Ray(const Vec3& _origin, const Vec3& _direction, double _time) 
  : origin(_origin), 
  direction(_direction),
	time(_time)
{
    direction.normalize();
}

void Ray::perturb(float roughness)
{
	Vec3 u, v;
	Vec3::createONB(direction, u, v);
	float psi1 = generateRandomFloat(0, 1);
	float psi2 = generateRandomFloat(0, 1);
	//std::cout << psi1 << psi2 << std::endl;
	direction = direction + (u * (psi1 - 0.5) + v * (psi2 - 0.5)) * roughness;
	direction.normalize();
}
