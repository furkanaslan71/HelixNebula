#include "mesh.h"

bool isDegenerate(const glm::vec3 indices[3]) {
	return false;
	const float eps = 1e-6f;
	const float epsSq = eps * eps;

	// 1. Check for Coincident Vertices (points too close together)
	// We calculate squared distance manually: dot(diff, diff)
	glm::vec3 d01 = indices[1] - indices[0];
	glm::vec3 d12 = indices[2] - indices[1];
	glm::vec3 d20 = indices[0] - indices[2];

	if (glm::dot(d01, d01) < epsSq ||
		glm::dot(d12, d12) < epsSq ||
		glm::dot(d20, d20) < epsSq)
	{
		return true;
	}

	// 2. Check for Collinearity (points forming a line, not a triangle)
	// A triangle is degenerate if its area is zero.
	// Area is proportional to the magnitude of the cross product of two edges.
	glm::vec3 crossProd = glm::cross(d01, indices[2] - indices[0]);

	// If squared magnitude of cross product is near zero, it's a line
	if (glm::dot(crossProd, crossProd) < epsSq)
	{
		return true;
	}

	return false;
}

template bool Mesh::hit<true>(const Ray&, Interval, HitRecord&) const;
template bool Mesh::hit<false>(const Ray&, Interval, HitRecord&) const;

Mesh::Mesh(int _id, const std::vector<Triangle_>& _faces,
	const std::vector<glm::vec3>& vertex_data,
	int vertex_offset, int texture_offset, bool is_smooth,
	const std::vector<glm::vec2>& tex_data,
	 std::optional<glm::vec3> _radiance)
	: id(_id), radiance(_radiance)
{
	if (is_smooth)
	{
		std::vector<std::vector<std::pair<glm::vec3, double>>> per_vertex_triangles; // pair<triangle_normal, area> for each vertex
		per_vertex_triangles.resize(vertex_data.size());
		for (const Triangle_& raw_triangle : _faces)
		{
			glm::vec3 v0 = glm::vec3(vertex_data[raw_triangle.v0_id + vertex_offset]);
			glm::vec3 v1 = glm::vec3(vertex_data[raw_triangle.v1_id + vertex_offset]);
			glm::vec3 v2 = glm::vec3(vertex_data[raw_triangle.v2_id + vertex_offset]);
			double area = Triangle::getAreaTriangle(v0, v1, v2);
			glm::vec3 edge1 = v1 - v0;
			glm::vec3 edge2 = v2 - v0;
			glm::vec3 face_normal = glm::normalize(glm::cross(edge1, edge2));
			per_vertex_triangles[raw_triangle.v0_id + vertex_offset].push_back(
				std::make_pair(face_normal, area)
			);
			per_vertex_triangles[raw_triangle.v1_id + vertex_offset].push_back(
				std::make_pair(face_normal, area)
			);
			per_vertex_triangles[raw_triangle.v2_id + vertex_offset].push_back(
				std::make_pair(face_normal, area)
			);
		}
		std::vector<glm::vec3> vertex_normals;
		for (const auto& v : per_vertex_triangles)
		{
			glm::vec3 normal(0.0, 0.0, 0.0);
			double total_area = 0.0;
			for (const auto& pair : v)
			{
				normal = normal + pair.first * (float) pair.second;
				total_area += pair.second;
			}
			if (total_area > 0.0)
			{
				normal = normal / (float) total_area;
				normal = glm::normalize(normal);
			}
			vertex_normals.push_back(normal);
		}
		for (const Triangle_& raw_triangle : _faces)
		{
			glm::vec3 indices[3] = {
				vertex_data[raw_triangle.v0_id + vertex_offset],
				vertex_data[raw_triangle.v1_id + vertex_offset],
				vertex_data[raw_triangle.v2_id + vertex_offset]
			};

			glm::vec3 per_vertex_normals[3] = {
				vertex_normals[raw_triangle.v0_id + vertex_offset],
				vertex_normals[raw_triangle.v1_id + vertex_offset],
				vertex_normals[raw_triangle.v2_id + vertex_offset]
			};

			glm::vec2 tex_coords[3];
			if (!tex_data.empty())
			{
				tex_coords[0] = tex_data[raw_triangle.v0_id + texture_offset];
				tex_coords[1] = tex_data[raw_triangle.v1_id + texture_offset];
				tex_coords[2] = tex_data[raw_triangle.v2_id + texture_offset];
			}
			else
			{
				tex_coords[0] = {0.f, 0.f};
				tex_coords[1] = {0.f, 0.f};
				tex_coords[2] = {0.f, 0.f};
			}

			if (isDegenerate(indices))
				continue;
			faces.push_back(Triangle(indices, per_vertex_normals, tex_coords));
		}
	}
	else
	{
		for (const Triangle_& raw_triangle : _faces)
		{
			glm::vec3 indices[3] = { vertex_data[raw_triangle.v0_id + vertex_offset],
				vertex_data[raw_triangle.v1_id + vertex_offset],
				vertex_data[raw_triangle.v2_id + vertex_offset] };

			glm::vec2 tex_coords[3];
			if (!tex_data.empty())
			{
				tex_coords[0] = tex_data[raw_triangle.v0_id + texture_offset];
				tex_coords[1] = tex_data[raw_triangle.v1_id + texture_offset];
				tex_coords[2] = tex_data[raw_triangle.v2_id + texture_offset];
			}
			else
			{
				tex_coords[0] = {0.f, 0.f};
				tex_coords[1] = {0.f, 0.f};
				tex_coords[2] = {0.f, 0.f};
			}
			if (isDegenerate(indices))
				continue;

			faces.push_back(Triangle(indices, tex_coords));
		}
	}
	bvh = std::make_shared<BVH<Triangle>>(faces);
	bounding_box = bvh->getAABB();
}

template <bool occlusion_only>
bool Mesh::hit(const Ray& ray, Interval ray_t, HitRecord& rec) const
{
	bool res = bvh->intersect<occlusion_only>(ray, ray_t, rec);
	if (res)
	{
		rec.radiance = radiance;
	}
	return res;
}

AABB Mesh::getAABB() const
{
	return bounding_box;
}