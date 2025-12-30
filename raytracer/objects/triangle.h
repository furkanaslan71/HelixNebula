#ifndef TRIANGLE_H
#define TRIANGLE_H

#include "core/hittable.h"
#include "parser/parser.hpp"
#include "core/aabb.h"
#include "glm_config.h"


class Triangle{
public:

	Triangle(glm::vec3 _indices[3])
		: indices{_indices[0], _indices[1], _indices[2]}
	{
		glm::vec3 min = indices[0];
		glm::vec3 max = indices[0];
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
		glm::vec3 vec1 = indices[1] - indices[0];
		glm::vec3 vec2 = indices[2] - indices[0];
		vec1 = glm::cross(vec1, vec2);
		vec1 = glm::normalize(vec1);
		this->normal = vec1;
	}

	Triangle(glm::vec3 _indices[3], glm::vec3 _per_vertex_normals[3])
		: 
		indices{ _indices[0], _indices[1], _indices[2] },
		smooth_shading(true), 
		per_vertex_normals{
			_per_vertex_normals[0], 
			_per_vertex_normals[1],
			_per_vertex_normals[2] 
		}
	{
		glm::vec3 min = indices[0];
		glm::vec3 max = indices[0];
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
		glm::vec3 vec1 = indices[1] - indices[0];
		glm::vec3 vec2 = indices[2] - indices[0];
		vec1 = glm::cross(vec1, vec2);
		vec1 = glm::normalize(vec1);
		this->normal = vec1;
	}


	template<bool occlusion_only>
	bool hit(const Ray& ray, Interval ray_t, HitRecord& rec) const
	{
		glm::vec3 c1 = indices[0] - indices[1];
		glm::vec3 c2 = indices[0] - indices[2];
		glm::vec3 c3 = ray.direction;
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
		if constexpr (occlusion_only) return true;

		if (beta + gamma <= 1 && beta + 0.00000001 >= 0 && gamma + 0.00000001 >= 0)
		{
			rec.t = t;
			rec.point = ray.origin + ray.direction * (float)t;
			if (this->smooth_shading)
			{
				glm::vec3 barycentric_coords = barycentricCoefficients(rec.point);
				rec.normal = per_vertex_normals[0] * barycentric_coords.x +
					per_vertex_normals[1] * barycentric_coords.y +
					per_vertex_normals[2] * barycentric_coords.z;
				rec.normal = glm::normalize(rec.normal);
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

	AABB getAABB() const { return bounding_box; }

	static inline double getAreaTriangle(glm::vec3 v1, glm::vec3 v2, glm::vec3 v3)
	{
		glm::vec3 edge1 = v2 - v1;
		glm::vec3 edge2 = v3 - v1;
		glm::vec3 cross_product = glm::cross(edge1, edge2);
		double area = 0.5 * glm::length(cross_product);
		return area;
	}

private:
	glm::vec3 normal;
	glm::vec3 indices[3];
	AABB bounding_box;
	bool smooth_shading = false;
	glm::vec3 per_vertex_normals[3];

	inline double det(const glm::vec3& c0, const glm::vec3& c1, const glm::vec3& c2) const
	{
		double temp1 = c0.x *
			(c1.y * c2.z - c1.z * c2.y);

		double temp2 = c1.x *
			(c0.y * c2.z - c0.z * c2.y);

		double temp3 = c2.x *
			(c0.y * c1.z - c0.z * c1.y);

		return temp1 - temp2 + temp3;
	}

	inline glm::vec3 barycentricCoefficients(const glm::vec3& point) const
	{
		glm::vec3 v0 = indices[1] - indices[0];
		glm::vec3 v1 = indices[2] - indices[0];
		glm::vec3 v2 = point - indices[0];
		double d00 = glm::dot(v0, v0);
		double d01 = glm::dot(v0, v1);
		double d11 = glm::dot(v1, v1);
		double d20 = glm::dot(v2, v0);
		double d21 = glm::dot(v2, v1);
		double denom = d00 * d11 - d01 * d01;
		double v = (d11 * d20 - d01 * d21) / denom;
		double w = (d00 * d21 - d01 * d20) / denom;
		double u = 1.0 - v - w;
		return glm::vec3(u, v, w);
	}

	
};

#endif // !TRIANGLE_H