#ifndef VECTOR_UTILITY_H
#define VECTOR_UTILITY_H
#include "external/glm/glm/glm.hpp"

static inline void createONB(const glm::vec3& r, glm::vec3& u, glm::vec3& v)
{
  glm::vec3 w = glm::normalize(r);

  glm::vec3 helper = (std::abs(w.x) > 0.9f)
    ? glm::vec3(0.0f, 1.0f, 0.0f)
    : glm::vec3(1.0f, 0.0f, 0.0f);

  u = glm::normalize(glm::cross(helper, w));

  v = glm::cross(w, u);
}

inline bool isZero(const glm::vec3& v)
{
  return (v.x == 0.0f && v.y == 0.0f && v.z == 0.0f);
}

#endif //!VECTOR_UTILITY_H