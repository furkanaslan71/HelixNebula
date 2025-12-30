#ifndef GEOMETRY_CONCEPTS_H
#define GEOMETRY_CONCEPTS_H
#include <concepts>

#include "core/ray.h"
#include "core/interval.h"
#include "core/hittable.h"
#include "core/aabb.h"

template <typename T>
concept HasHit = requires(const T & obj, const Ray & ray, Interval ray_t, HitRecord & rec)
{
  { obj.template hit<false>(ray, ray_t, rec) } -> std::same_as<bool>;
  { obj.template hit<true>(ray, ray_t, rec) } -> std::same_as<bool>;
};

template <typename T>
concept HasGetAABB = requires(const T & obj)
{
  { obj.getAABB() } -> std::same_as<AABB>;
};

template <typename T>
concept GeometryConcept = HasHit<T> && HasGetAABB<T>;

#endif // !GEOMETRY_CONCEPTS_H
