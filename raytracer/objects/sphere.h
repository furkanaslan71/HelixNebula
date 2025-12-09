#ifndef SPHERE_H
#define SPHERE_H

#include "core/hittable.h"
#include "parser/parser.hpp"
#include "core/aabb.h"


class Sphere{
public:
	glm::vec3 center;
	double radius;

	Sphere(glm::vec3 _center, double _radius)
		: 
		center(_center), 
		radius(_radius)
	{
		bounding_box = AABB(center - glm::vec3(radius, radius, radius),
														center + glm::vec3(radius, radius, radius));
	}

	bool hit(const Ray& ray, Interval ray_t, HitRecord& rec) const
	{

		glm::vec3 oc = ray.origin - center;
		double a = glm::dot(ray.direction, ray.direction);
		double b = glm::dot(oc, ray.direction);
		double c = glm::dot(oc, oc) - radius * radius;
		double discriminant = b * b - a * c;

		// Intersection occurs
		if (discriminant >= 0)
		{
			double t1 = (-b - sqrt(discriminant)) / a;
			double t2 = (-b + sqrt(discriminant)) / a;

			if (t1 >= 0 || t2 >= 0)
			{
				double t = (t1 >= 0) ? t1 : t2;

				glm::vec3 hitPoint = ray.origin + ray.direction * (float)t;
				glm::vec3 normal = glm::normalize((hitPoint - center));

				rec.t = t;
				rec.point = hitPoint;
				rec.normal = normal;
				rec.set_front_face(ray);
				return true;
			}
		}
		return false;

	};
	AABB getAABB() const
	{
		return bounding_box;
	}

private:
	AABB bounding_box;
};

#endif // SPHERE_H