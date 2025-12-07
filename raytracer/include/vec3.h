#ifndef VEC3_H
#define VEC3_H
#include <iostream>
#include <cmath>
#include "parser.hpp"
#include <../external/glm/glm/glm.hpp>
#include <../external/glm/glm/gtc/matrix_transform.hpp>

class Vec3 {
public:
    double x, y, z;

    Vec3() : x(0), y(0), z(0) {}
    Vec3(double _x, double _y, double _z) : x(_x), y(_y), z(_z) {}
		Vec3(double val) : x(val), y(val), z(val) {}
		Vec3(const Vec3f_& vec) : x(static_cast<double>(vec.x)), y(static_cast<double>(vec.y)), z(static_cast<double>(vec.z)) {}
    Vec3(const glm::vec3 glm_vec)
    {
      x = static_cast<double>(glm_vec.x);
      y = static_cast<double>(glm_vec.y);
      z = static_cast<double>(glm_vec.z);
    }
    Vec3(const glm::vec4 glm_vec)
    {
      x = static_cast<double>(glm_vec.x);
      y = static_cast<double>(glm_vec.y);
      z = static_cast<double>(glm_vec.z);
    }

    inline glm::vec3 toGlm3()
    {
      return glm::vec3(
        static_cast<float>(x),
        static_cast<float>(y),
        static_cast<float>(z)
      );
    }

    inline glm::vec4 toGlm4()
    {
      return glm::vec4(
        static_cast<float>(x),
        static_cast<float>(y),
        static_cast<float>(z),
        1.0f
      );
    }
    
    // Mutating
    inline Vec3 normalize() {
        double len = length();
        if (len < 1e-8)
          len = 1e-8;
        x /= len;
        y /= len;
        z /= len;
        return *this;
    }
    inline Vec3 reverse()
    {
      x *= -1;
      y *= -1;
      z *= -1;
      return *this;
    }
    
    // Non-mutating
    inline Vec3 inverse() const
    {
      return Vec3(x * -1, y * -1, x * -1);
    }
    inline Vec3 normalized() const
    {
      double len = length();
      if (len < 1e-8)
        len = 1e-8;
      return Vec3(x / len, y / len, z / len);
    }

    inline double distance(const Vec3& _other) const
    {
      return std::sqrt(
        (x - _other.x) * (x - _other.x) +
        (y - _other.y) * (y - _other.y) +
				(z - _other.z) * (z - _other.z));
    }

    inline Vec3 cross(const Vec3& _other) const{
        return Vec3(
            y * _other.z - z * _other.y,
            z * _other.x - x * _other.z,
            x * _other.y - y * _other.x);
    }

    inline double dot(const Vec3& _other) const 
    {
        return x * _other.x + y * _other.y + z * _other.z;
    }

    inline double length() const 
    {
        return sqrt(x * x + y * y + z * z);
    }

    // Non-mutating operators
    inline Vec3 operator+(const Vec3& other) const
    {
      return Vec3(x + other.x, y + other.y, z + other.z);
    }

    inline Vec3 operator-(const Vec3& other) const
    {
      return Vec3(x - other.x, y - other.y, z - other.z);
    }

    inline Vec3 operator*(const Vec3& other) const
    {
      return Vec3(x * other.x, y * other.y, z * other.z);
    }

    inline Vec3 operator*(const double& scalar) const
    {
      return Vec3(x * scalar, y * scalar, z * scalar);
    }

    inline Vec3 operator/(const Vec3& other) const
    {
      return Vec3(x / other.x, y / other.y, z / other.z);
    }

    inline Vec3 operator/(const double& scalar) const
    {
      return Vec3(x / scalar, y / scalar, z / scalar);
    }

    inline double operator[](int i) const
    {
      if (i == 0) return x;
      if (i == 1) return y;
      if (i == 2) return z;
      throw std::out_of_range("vec3, axis index out of range");
    }

    // Mutating operators
    inline Vec3 operator+=(const Vec3& other)
    {
      x += other.x;
      y += other.y;
      z += other.z;
      return *this;
    }

    inline Vec3 operator-=(const Vec3& other)
    {
      x -= other.x;
      y -= other.y;
      z -= other.z;
      return *this;
    }

    inline Vec3 operator*=(const Vec3& other)
    {
      x *= other.x;
      y *= other.y;
      z *= other.z;
      return *this;
    }

    inline Vec3 operator*=(const double& scalar)
    {
      x *= scalar;
      y *= scalar;
      z *= scalar;
      return *this;
    }

    inline Vec3 operator/=(const Vec3& other)
    {
      x /= other.x;
      y /= other.y;
      z /= other.z;
      return *this;
    }

    inline Vec3 operator/=(const double& scalar)
    {
      x *= scalar;
      y *= scalar;
      z *= scalar;
      return *this;
    }

    inline bool isZero() const
    {
      return ((x == 0) && (y == 0) && (z == 0));
    }

    static inline void createONB(const  glm::vec3& r, glm::vec3& u, glm::vec3& v)
    {
       glm::vec3 w = glm::normalize(r);

       glm::vec3 helper = (std::abs(w.x) > 0.9f)
        ? glm::vec3(0.0f, 1.0f, 0.0f)
        : glm::vec3(1.0f, 0.0f, 0.0f);

      u = glm::normalize(glm::cross(helper, w));

      v = glm::cross(w, u);
    }

};

#endif // !VEC3_H