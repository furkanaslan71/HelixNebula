#ifndef TRIANGLE_H
#define TRIANGLE_H

#include "../include/hittable.h"
#include "../include/parser.hpp"
#include "../include/aabb.h"


class Triangle : public Hittable {
public:

	Triangle(Vec3 _indices[3])
		: indices{_indices[0], _indices[1], _indices[2]}
	{
		Vec3 min = indices[0];
		Vec3 max = indices[0];
		for (int i = 1; i < 3; i++)
		{
			min.x = fmin(indices[i].x, min.x);
			max.x = fmax(indices[i].x, max.x);
			min.y = fmin(indices[i].y, min.y);
			max.y = fmax(indices[i].y, max.y);
			min.z = fmin(indices[i].z, min.z);
			max.z = fmax(indices[i].z, max.z);
		}
		bounding_box = AABB(min, max);
		Vec3 vec1 = indices[1] - indices[0];
		Vec3 vec2 = indices[2] - indices[0];
		vec1 = vec1.cross(vec2);
		vec1.normalize();
		this->normal = vec1;
	}

	Triangle(Vec3 _indices[3], Vec3 _per_vertex_normals[3])
		: 
		indices{ _indices[0], _indices[1], _indices[2] },
		smooth_shading(true), 
		per_vertex_normals{
			_per_vertex_normals[0], 
			_per_vertex_normals[1],
			_per_vertex_normals[2] 
		}
	{
		Vec3 min = indices[0];
		Vec3 max = indices[0];
		for (int i = 1; i < 3; i++)
		{
			min.x = fmin(indices[i].x, min.x);
			max.x = fmax(indices[i].x, max.x);
			min.y = fmin(indices[i].y, min.y);
			max.y = fmax(indices[i].y, max.y);
			min.z = fmin(indices[i].z, min.z);
			max.z = fmax(indices[i].z, max.z);
		}
		bounding_box = AABB(min, max);
		Vec3 vec1 = indices[1] - indices[0];
		Vec3 vec2 = indices[2] - indices[0];
		vec1 = vec1.cross(vec2);
		vec1.normalize();
		this->normal = vec1;
	}



	bool hit(const Ray& ray, Interval ray_t, HitRecord& rec) const override
	{
		Vec3 c1 = indices[0] - indices[1];
		Vec3 c2 = indices[0] - indices[2];
		Vec3 c3 = ray.direction;
		double detA = det(c1, c2, c3);
		if (detA == 0) return false;

		c1 = indices[0] - ray.origin;
		double beta = det(c1, c2, c3) / detA;

		c2 = c1;
		c1 = indices[0] - indices[1];
		double gamma = det(c1, c2, c3) / detA;

		c3 = c2;
		c2 = indices[0] - indices[2];
		double t = det(c1, c2, c3) / detA;

		if (t < ray_t.min + 0.00000001 || 0.00000001 + t > ray_t.max) return false;

		if (beta + gamma <= 1 && beta + 0.00000001 >= 0 && gamma + 0.00000001 >= 0)
		{
			rec.t = t;
			rec.point = ray.origin + ray.direction * t;
			if (this->smooth_shading)
			{
				Vec3 barycentric_coords = barycentricCoefficients(rec.point);
				rec.normal = per_vertex_normals[0] * barycentric_coords.x +
					per_vertex_normals[1] * barycentric_coords.y +
					per_vertex_normals[2] * barycentric_coords.z;
				rec.normal.normalize();
			}
			else
			{
				rec.normal = this->normal;
			}
			rec.set_front_face(ray);
			return true;
		}
		return false;
		
	}

	AABB getAABB() const override { return bounding_box; }

	static inline double getAreaTriangle(Vec3 v1, Vec3 v2, Vec3 v3)
	{
		Vec3 edge1 = v2 - v1;
		Vec3 edge2 = v3 - v1;
		Vec3 cross_product = edge1.cross(edge2);
		double area = 0.5 * cross_product.length();
		return area;
	}

private:
	Vec3 normal;
	Vec3 indices[3];
	AABB bounding_box;
	bool smooth_shading = false;
	Vec3 per_vertex_normals[3];

	inline double det(const Vec3& c0, const Vec3& c1, const Vec3& c2) const
	{
		double temp1 = c0.x *
			(c1.y * c2.z - c1.z * c2.y);

		double temp2 = c1.x *
			(c0.y * c2.z - c0.z * c2.y);

		double temp3 = c2.x *
			(c0.y * c1.z - c0.z * c1.y);

		return temp1 - temp2 + temp3;
	}

	inline Vec3 barycentricCoefficients(const Vec3& point) const
	{
		Vec3 v0 = indices[1] - indices[0];
		Vec3 v1 = indices[2] - indices[0];
		Vec3 v2 = point - indices[0];
		double d00 = v0.dot(v0);
		double d01 = v0.dot(v1);
		double d11 = v1.dot(v1);
		double d20 = v2.dot(v0);
		double d21 = v2.dot(v1);
		double denom = d00 * d11 - d01 * d01;
		double v = (d11 * d20 - d01 * d21) / denom;
		double w = (d00 * d21 - d01 * d20) / denom;
		double u = 1.0 - v - w;
		return Vec3(u, v, w);
	}

	
};

#endif // !TRIANGLE_H