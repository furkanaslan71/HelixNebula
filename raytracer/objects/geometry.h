#ifndef GEOMETRY_H
#define GEOMETRY_H
#include <variant>
#include <concepts>

#include "geometry_concepts.h"
#include "mesh.h"
#include "triangle.h"
#include "sphere.h"


template <typename Variant>
struct validate_variant; // primary template (unspecialized)

template <typename... Ts>
struct validate_variant<std::variant<Ts...>> { // specialization for std::variant
  static_assert((GeometryConcept<Ts> && ...),
                "All types in Geometry::Variant must satisfy GeometryConcept");
};


class Geometry {
public:
	using Variant = std::variant<Triangle, Sphere, Mesh>;

  static validate_variant<Variant> _check;

  template <GeometryConcept T, typename... Args>
  Geometry(std::in_place_type_t<T>, Args&&... args)
    : data(std::in_place_type<T>, std::forward<Args>(args)...)
  {
  }

  /*template <typename T, typename... Args>
  T& emplace(Args&&... args)
  {
    return data.emplace<T>(std::forward<Args>(args)...);
  }*/

  bool hit(const Ray& ray, Interval ray_t, HitRecord& rec) const
  {
    return std::visit([&](auto const& obj) -> bool {
      return obj.hit(ray, ray_t, rec);
                      }, data);
  }

  AABB getAABB() const
  {
    return std::visit([&](auto const& obj) -> AABB {
      return obj.getAABB();
                      }, data);
  }

private:
  Variant data;

};

template <GeometryConcept T>
class GeometryTemplated {
public:
  template <typename... Args>
  explicit GeometryTemplated(Args&&... args)
    : data(std::forward<Args>(args)...)
  { }

  bool hit(const Ray& ray, Interval ray_t, HitRecord& rec) const
  {
    return data.hit(ray, ray_t, rec);
  }

  AABB getAABB() const
  {
    return data.getAABB();
  }

private:
  T data;
};



#endif //!GEOMETRY_H