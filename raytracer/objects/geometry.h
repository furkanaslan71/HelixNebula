#ifndef GEOMETRY_H
#define GEOMETRY_H
#include <variant>

#include "mesh.h"
#include "triangle.h"
#include "sphere.h"
#include "tlas_box.h"

template <typename T>
concept HasHit = requires(const T & obj, const Ray & ray, Interval ray_t, HitRecord & rec)
{
  { obj.hit(ray, ray_t, rec) } -> std::same_as<bool>;
};

template <typename T>
concept HasGetAABB = requires(const T & obj)
{
  { obj.getAABB() } -> std::same_as<AABB>;
};

template <typename T>
concept GeometryConcept = HasHit<T> && HasGetAABB<T>;

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




#endif //!GEOMETRY_H