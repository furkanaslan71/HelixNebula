#ifndef SPHERE_H
#define SPHERE_H
#include <numbers>
#include "core/hittable.h"
#include "parser/parser.hpp"
#include "core/aabb.h"


class Sphere{
public:
	glm::vec3 center;
	double radius;
	bool texture;

	Sphere(glm::vec3 _center, double _radius, bool _texture = false)
		: 
		center(_center), 
		radius(_radius),
		texture(_texture)
	{
		bounding_box = AABB(center - glm::vec3(radius, radius, radius),
														center + glm::vec3(radius, radius, radius));
	}

	void getSphereUV(const glm::vec3& hitPoint, float& u, float& v) const
	{
		float pi = glm::pi<float>();
		auto P = hitPoint - center;
		auto theta = std::acos(P.y / radius	);
		auto phi = std::atan2(P.z, P.x);
		
		u = (-phi + pi) / (2.0f * pi);
		v = theta / pi;
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