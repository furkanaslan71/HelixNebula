#ifndef GEOMETRY_CONCEPTS_H
#define GEOMETRY_CONCEPTS_H
#include <concepts>

#include "../include/ray.h"
#include "../include/interval.h"
#include "../include/hittable.h"
#include "../include/aabb.h"

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

#endif // !GEOMETRY_CONCEPTS_H
