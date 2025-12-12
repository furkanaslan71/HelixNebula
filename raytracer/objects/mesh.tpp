
template<Shading mode, TextureLookup tex>
Mesh<mode, tex>::Mesh(int _id, const std::vector<Triangle_>& _faces,
	const std::vector<glm::vec3>& vertex_data, const std::vector<glm::vec2> uv_data)
	: id(_id)
{
	if constexpr (mode == Shading::Smooth)
	{
		std::vector<std::vector<std::pair<glm::vec3, double>>> per_vertex_triangles; // pair<triangle_normal, area> for each vertex
		per_vertex_triangles.resize(vertex_data.size());
		for (const Triangle_& raw_triangle : _faces)
		{
			glm::vec3 v0 = glm::vec3(vertex_data[raw_triangle.v0_id]);
			glm::vec3 v1 = glm::vec3(vertex_data[raw_triangle.v1_id]);
			glm::vec3 v2 = glm::vec3(vertex_data[raw_triangle.v2_id]);
			double area = Triangle::getAreaTriangle(v0, v1, v2);
			glm::vec3 edge1 = v1 - v0;
			glm::vec3 edge2 = v2 - v0;
			glm::vec3 face_normal = glm::normalize(glm::cross(edge1, edge2));
			per_vertex_triangles[raw_triangle.v0_id].push_back(
				std::make_pair(face_normal, area)
			);
			per_vertex_triangles[raw_triangle.v1_id].push_back(
				std::make_pair(face_normal, area)
			);
			per_vertex_triangles[raw_triangle.v2_id].push_back(
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
			glm::vec3 indices[3] = { vertex_data[raw_triangle.v0_id],
				vertex_data[raw_triangle.v1_id],
				vertex_data[raw_triangle.v2_id] };

			PerVertexNormals per_vertex_normals = {
				vertex_normals[raw_triangle.v0_id],
				vertex_normals[raw_triangle.v1_id],
				vertex_normals[raw_triangle.v2_id] };

			if constexpr (tex == TextureLookup::Textured)
			{
				auto tex_coords = TexCoords(raw_triangle, uv_data);
				faces.push_back(TriangleNew<mode, tex>(indices[0], indices[1], indices[2], per_vertex_normals, tex_coords));
			}
			else
			{
				faces.push_back(TriangleNew<mode, tex>(indices[0], indices[1], indices[2], per_vertex_normals));
			}
		}
	}
	else
	{
		for (const Triangle_& raw_triangle : _faces)
		{
			glm::vec3 indices[3] = { vertex_data[raw_triangle.v0_id],
				vertex_data[raw_triangle.v1_id],
				vertex_data[raw_triangle.v2_id] };

			if constexpr (tex == TextureLookup::Textured)
			{
				auto tex_coords = TexCoords(raw_triangle, uv_data);
				faces.push_back(TriangleNew<mode, tex>(indices[0], indices[1], indices[2], tex_coords));
			}
			else
			{
				faces.push_back(TriangleNew<mode, tex>(indices[0], indices[1], indices[2]));
			}
		}
	}
	bvh = std::make_shared<BVH<TriangleNew<mode, tex>>>(faces);
	bounding_box = bvh->getAABB();
}

template<Shading mode, TextureLookup tex>
bool Mesh<mode, tex>::hit(const Ray& ray, Interval ray_t, HitRecord& rec) const
{
	return bvh->intersect(ray, ray_t, rec);
}
template<Shading mode, TextureLookup tex>
AABB Mesh<mode, tex>::getAABB() const
{
	return bounding_box;
}