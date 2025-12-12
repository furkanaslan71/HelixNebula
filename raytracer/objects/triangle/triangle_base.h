#ifndef TRIANGLE_BASE_H
#define TRIANGLE_BASE_H
#include "glm_config.h"
#include "core/hittable.h"
#include "core/aabb.h"
#include "core/ray.h"
#include "core/interval.h"

enum class Shading : bool {
  Flat = false,
  Smooth = true
};

enum class TextureLookup : bool {
  NoTexture = false,
  Textured = true
};

struct TexCoords {
  glm::vec2 uvs[3];
  TexCoords()
  {
    uvs[0] = glm::vec2{-1.f};
    uvs[1] = glm::vec2{ -1.f };
    uvs[2] = glm::vec2{ -1.f };
  }
  TexCoords(const Triangle_& tri, const std::vector<glm::vec2>& uv_data)
  {
    uvs[0] = uv_data[tri.v0_id];
    uvs[1] = uv_data[tri.v1_id];
    uvs[2] = uv_data[tri.v2_id];
  };
};

struct BarycentricModule {
  double d00, d01, d11, denom;

  BarycentricModule()
  {
    d00 = 0;
    d01 = 0;
    d11 = 0;
    denom = 0;
  }

  BarycentricModule(const glm::vec3& e1, const glm::vec3& e2)
  {
    glm::vec3 i0 = e1;
    glm::vec3 i1 = e2;
    d00 = glm::dot(i0, i0);
    d01 = glm::dot(i0, i1);
    d11 = glm::dot(i1, i1);
    denom = d00 * d11 - d01 * d01;
  }

  inline glm::vec3 getBarycentricCoefficients(
    const glm::vec3& point,
    const glm::vec3& v0,
    const glm::vec3& e1,
    const glm::vec3& e2) const
  {
    glm::vec3 a0 = e1;
    glm::vec3 a1 = e2;
    glm::vec3 a2 = point - v0;
    float d20 = glm::dot(a2, a0);
    float d21 = glm::dot(a2, a1);
    float v = (d11 * d20 - d01 * d21) / denom;
    float w = (d00 * d21 - d01 * d20) / denom;
    float u = 1.0 - v - w;
    return glm::vec3(u, v, w);
  }
};

struct PerVertexNormals {
  glm::vec3 n1, n2, n3;
};

template<Shading mode, TextureLookup tex>
class TriangleNew {
public:
  template<typename... Extra>
  TriangleNew(const glm::vec3& a,
               const glm::vec3& b,
               const glm::vec3& c,
              Extra&&... extra)
    : v0(a), v1(b), v2(c)
  {
    e1 = v1 - v0;
    e2 = v2 - v0;
    normal = glm::normalize(glm::cross(e1, e2));

    glm::vec3 min;
    glm::vec3 max;
    min.x = fmin(fmin(v0.x, v1.x), v2.x);
    min.y = fmin(fmin(v0.y, v1.y), v2.y);
    min.z = fmin(fmin(v0.z, v1.z), v2.z);
    max.x = fmax(fmax(v0.x, v1.x), v2.x);
    max.y = fmax(fmax(v0.y, v1.y), v2.y);
    max.z = fmax(fmax(v0.z, v1.z), v2.z);
    bounding_box = AABB(min, max);

    initExtra(std::forward<Extra>(extra)...);
  }

  bool hit(const Ray& ray, Interval ray_t, HitRecord& rec) const
  {
    glm::vec3 c1 = -e1;
    glm::vec3 c2 = -e2;
    glm::vec3 c3 = ray.direction;
    double detA = det(c1, c2, c3);
    if (fabs(detA) < 1e-12) return false;

    c1 = v0 - ray.origin;
    double beta = det(c1, c2, c3) / detA;

    c2 = c1;
    c1 = -e1;
    double gamma = det(c1, c2, c3) / detA;

    c3 = c2;
    c2 = -e2;
    double t = det(c1, c2, c3) / detA;

    constexpr double EPS = 1e-8;

    if (t < ray_t.min + EPS || EPS + t > ray_t.max) return false;

    if (beta + gamma <= 1 && beta + EPS >= 0 && gamma + EPS >= 0)
    {
      rec.t = t;
      rec.point = ray.origin + ray.direction * (float)t;
      if constexpr (mode == Shading::Flat && tex == TextureLookup::NoTexture)
      {
        rec.normal = normal;
      }
      else if constexpr (mode == Shading::Smooth && tex == TextureLookup::NoTexture)
      {
        glm::vec3 barycentric = bary.getBarycentricCoefficients(
          rec.point, v0, e1, e2);
        rec.normal = per_vertex_normals.n1 * barycentric.x +
          per_vertex_normals.n2 * barycentric.y +
          per_vertex_normals.n3 * barycentric.z;
        rec.normal = glm::normalize(rec.normal);
      }
      else if constexpr (mode == Shading::Flat && tex == TextureLookup::Textured)
      {
        glm::vec3 barycentric = bary.getBarycentricCoefficients(
          rec.point, v0, e1, e2);
        rec.uv = tex_coords.uvs[0] * barycentric.x +
          tex_coords.uvs[1] * barycentric.y +
          tex_coords.uvs[2] * barycentric.z;
      }
      else if constexpr (mode == Shading::Smooth && tex == TextureLookup::Textured)
      {
        glm::vec3 barycentric = bary.getBarycentricCoefficients(
          rec.point, v0, e1, e2);
        rec.normal = per_vertex_normals.n1 * barycentric.x +
          per_vertex_normals.n2 * barycentric.y +
          per_vertex_normals.n3 * barycentric.z;
        rec.normal = glm::normalize(rec.normal);
        rec.uv = tex_coords.uvs[0] * barycentric.x +
          tex_coords.uvs[1] * barycentric.y +
          tex_coords.uvs[2] * barycentric.z;
      }
      rec.set_front_face(ray);
      return true;      
    }
    return false;

  }

  AABB getAABB() const { return bounding_box; }

  AABB bounding_box;
  glm::vec3 v0, v1, v2;
  glm::vec3 e1, e2;
  glm::vec3 normal;

  std::conditional_t<mode == Shading::Smooth,
    PerVertexNormals, std::monostate> per_vertex_normals;

  std::conditional_t<tex == TextureLookup::Textured,
    TexCoords, std::monostate> tex_coords;

  static constexpr bool needs_bary =
    (mode == Shading::Smooth) || (tex == TextureLookup::Textured);

  std::conditional_t<needs_bary,
    BarycentricModule, std::monostate> bary;

private:
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

  template<typename... Extra>
  void initExtra(Extra&&... extra)
  {
    if constexpr (mode == Shading::Smooth && tex == TextureLookup::Textured)
    {
      static_assert(sizeof...(Extra) == 2,
                    "Smooth + Textured triangles need (PerVertexNormals, TexCoords)");
      auto tuple = std::forward_as_tuple(extra...);
      per_vertex_normals = std::get<0>(tuple);
      tex_coords = std::get<1>(tuple);
      bary = BarycentricModule(e1, e2);
    }
    else if constexpr (mode == Shading::Smooth)
    {
      static_assert(sizeof...(Extra) == 1,
                    "Smooth triangles require PerVertexNormals");
      per_vertex_normals = std::get<0>(std::forward_as_tuple(extra...));
      bary = BarycentricModule(e1, e2);
    }
    else if constexpr (tex == TextureLookup::Textured)
    {
      static_assert(sizeof...(Extra) == 1,
                    "Textured triangles require TexCoords");
      tex_coords = std::get<0>(std::forward_as_tuple(extra...));
      bary = BarycentricModule(e1, e2);
    }
    else
    {
      static_assert(sizeof...(Extra) == 0,
                    "Flat / NoTexture triangles need no extra args");
    }
  }
};



#endif // !TRIANGLE_BASE_H