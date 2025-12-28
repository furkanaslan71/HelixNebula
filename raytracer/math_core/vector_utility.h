#ifndef VECTOR_UTILITY_H
#define VECTOR_UTILITY_H
#include <xmmintrin.h>
#include <iostream>

#include "glm_config.h"

/*
alignas(16) static const __m128 HALF = _mm_set1_ps(0.5f);
alignas(16) static const __m128 THREE_HALFS = _mm_set1_ps(1.5f);

static const __m256 half = _mm256_set1_ps(0.5f);
static const __m256 three_halfs = _mm256_set1_ps(1.5f);
*/

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

/*
inline glm::vec3 fastNormalizeSSE_NR(const glm::vec3& v)
{
  __m128 x = _mm_set_ps(0, v.z, v.y, v.x);
  __m128 dot = _mm_dp_ps(x, x, 0x7F);
  __m128 y = _mm_rsqrt_ps(dot);

  // Newton–Raphson
  __m128 y2 = _mm_mul_ps(y, y);
  __m128 x_y2 = _mm_mul_ps(dot, y2);
  y = _mm_mul_ps(y, _mm_sub_ps(THREE_HALFS, _mm_mul_ps(HALF, x_y2)));

  __m128 out = _mm_mul_ps(x, y);

  // Extract components WITHOUT array
  float X = _mm_cvtss_f32(out);
  float Y = _mm_cvtss_f32(_mm_shuffle_ps(out, out, 0x55));
  float Z = _mm_cvtss_f32(_mm_shuffle_ps(out, out, 0xAA));

  return { X, Y, Z };
  float result[4];
  _mm_storeu_ps(result, out); // Store the entire __m128 register to memory
  return { result[0], result[1], result[2] };
}

inline glm::mat4 fastMatrixMultiply(const glm::mat4& A, const glm::mat4& B)
{
  glm::mat4 C;

  // Load A columns (no transpose!)
  __m128 a0 = _mm_loadu_ps(&A[0][0]);
  __m128 a1 = _mm_loadu_ps(&A[1][0]);
  __m128 a2 = _mm_loadu_ps(&A[2][0]);
  __m128 a3 = _mm_loadu_ps(&A[3][0]);

  for (int col = 0; col < 4; ++col)
  {
    const float* bc = &B[col][0];

    __m128 b0 = _mm_broadcast_ss(bc + 0);
    __m128 b1 = _mm_broadcast_ss(bc + 1);
    __m128 b2 = _mm_broadcast_ss(bc + 2);
    __m128 b3 = _mm_broadcast_ss(bc + 3);

    // C[col] = a0*b0 + a1*b1 + a2*b2 + a3*b3
    __m128 c = _mm_mul_ps(a0, b0);
    c = _mm_fmadd_ps(a1, b1, c);
    c = _mm_fmadd_ps(a2, b2, c);
    c = _mm_fmadd_ps(a3, b3, c);

    _mm_storeu_ps(&C[col][0], c);
  }

  return C;
}


inline void fastMat4Vec4_2_AVX2(
  const glm::mat4& M,
  const glm::vec4& v0,
  const glm::vec4& v1,
  glm::vec4& out0,
  glm::vec4& out1)
{
  // --- 1. Load Matrix Columns (M[0] to M[3]) ---
  // GLM stores matrices in column-major order. M[c] is the c-th column.
  // We load them as 128-bit vectors.
  const __m128 m0 = _mm_load_ps(&M[0][0]); // Column 0 (x_x, y_x, z_x, w_x)
  const __m128 m1 = _mm_load_ps(&M[1][0]); // Column 1 (x_y, y_y, z_y, w_y)
  const __m128 m2 = _mm_load_ps(&M[2][0]); // Column 2 (x_z, y_z, z_z, w_z)
  const __m128 m3 = _mm_load_ps(&M[3][0]); // Column 3 (x_w, y_w, z_w, w_w)

  // --- 2. Interleave Input Vectors (v0 and v1) into a single 256-bit register ---
  // We load v0 and v1 into a 256-bit register, resulting in:
  // v_combined = [v0.x, v0.y, v0.z, v0.w | v1.x, v1.y, v1.z, v1.w]
  const __m256 v_combined = _mm256_set_m128(
    _mm_load_ps(&v1.x), // High 128 bits
    _mm_load_ps(&v0.x)  // Low 128 bits
  );

  // --- 3. Extract and Broadcast Components from v_combined ---
  // We need 4 broadcasted vectors for the X, Y, Z, W components of v0/v1 pair.

  // v_x_bcast = [v0.x, v0.x, v0.x, v0.x | v1.x, v1.x, v1.x, v1.x]
  // _MM_SHUFFLE(3, 3, 3, 3) is (7, 6, 5, 4) in the high 128 bits, (3, 2, 1, 0) in the low 128 bits.
  // _MM_SHUFFLE(0, 0, 0, 0) is picking index 0 (v0.x) from low and (v1.x) from high.
  const __m256 v_x_bcast = _mm256_shuffle_ps(v_combined, v_combined, _MM_SHUFFLE(0, 0, 0, 0));

  // v_y_bcast = [v0.y, v0.y, v0.y, v0.y | v1.y, v1.y, v1.y, v1.y]
  const __m256 v_y_bcast = _mm256_shuffle_ps(v_combined, v_combined, _MM_SHUFFLE(1, 1, 1, 1));

  // v_z_bcast = [v0.z, v0.z, v0.z, v0.z | v1.z, v1.z, v1.z, v1.z]
  const __m256 v_z_bcast = _mm256_shuffle_ps(v_combined, v_combined, _MM_SHUFFLE(2, 2, 2, 2));

  // v_w_bcast = [v0.w, v0.w, v0.w, v0.w | v1.w, v1.w, v1.w, v1.w]
  const __m256 v_w_bcast = _mm256_shuffle_ps(v_combined, v_combined, _MM_SHUFFLE(3, 3, 3, 3));


  // --- 4. Matrix-Vector Multiplication (Simultaneous for v0 and v1) ---
  // The matrix multiplication is: out = M[0]*v.x + M[1]*v.y + M[2]*v.z + M[3]*v.w
  // We combine the 128-bit matrix columns into 256-bit vectors where:
  // m_256_i = [M[i].x, M[i].y, M[i].z, M[i].w | M[i].x, M[i].y, M[i].z, M[i].w]

  // Create 256-bit matrix column vectors by duplicating the 128-bit column
  const __m256 m0_256 = _mm256_insertf128_ps(_mm256_castps128_ps256(m0), m0, 1);
  const __m256 m1_256 = _mm256_insertf128_ps(_mm256_castps128_ps256(m1), m1, 1);
  const __m256 m2_256 = _mm256_insertf128_ps(_mm256_castps128_ps256(m2), m2, 1);
  const __m256 m3_256 = _mm256_insertf128_ps(_mm256_castps128_ps256(m3), m3, 1);


  // Calculate term 1: M[0]*v.x
  // This calculates: [M[0]*v0.x | M[0]*v1.x]
  __m256 res = _mm256_mul_ps(m0_256, v_x_bcast);

  // Calculate term 2: M[1]*v.y and add it to the result
  // res += M[1]*v.y
  res = _mm256_fmadd_ps(m1_256, v_y_bcast, res); // Fused Multiply-Add

  // Calculate term 3: M[2]*v.z and add it to the result
  // res += M[2]*v.z
  res = _mm256_fmadd_ps(m2_256, v_z_bcast, res);

  // Calculate term 4: M[3]*v.w and add it to the result
  // res += M[3]*v.w
  res = _mm256_fmadd_ps(m3_256, v_w_bcast, res);


  // --- 5. Store Results ---
  // The result 'res' is [out0.x, out0.y, out0.z, out0.w | out1.x, out1.y, out1.z, out1.w]

  // Extract the low 128 bits (out0) and store
  _mm_store_ps(&out0.x, _mm256_castps256_ps128(res));

  // Extract the high 128 bits (out1) and store
  _mm_store_ps(&out1.x, _mm256_extractf128_ps(res, 1));
}

inline void printMat4(const glm::mat4& M)
{
  for (int r = 0; r < 4; ++r)
  {
    std::cout << "[ ";
    for (int c = 0; c < 4; ++c)
      std::cout << M[c][r] << " ";
    std::cout << "]\n";
  }
}



#ifdef __AVX2__
inline glm::vec3 fastNormalizeAVX(const glm::vec3& v)
{
  // Load x, y, z, 0 into a 256-bit register
  __m256 x = _mm256_set_ps(0, 0, 0, 0, 0.0f, v.z, v.y, v.x);

  // dot = x*x + y*y + z*z  (manual dot product, faster than _mm_dp_ps)
  __m256 sq = _mm256_mul_ps(x, x);

  // Horizontal add: sum first 3 elements (x*x + y*y + z*z)
  __m128 lo = _mm256_castps256_ps128(sq);                       // [x² y² z² 0]
  __m128 sh = _mm_movehdup_ps(lo);                              // [y² y² 0 0]
  __m128 sum1 = _mm_add_ps(lo, sh);                             // x²+y² | y²+z² | z²+0 | 0
  __m128 sh2 = _mm_movehl_ps(sh, sum1);                         // [z²+0 | 0 | ? | ?]
  __m128 dot = _mm_add_ss(sum1, sh2);                           // dot in lowest lane

  // Replicate dot across 256-bit lane (needed for AVX math)
  __m256 dot256 = _mm256_broadcastss_ps(dot);

  // initial Y = 1 / sqrt(dot)
  __m256 y = _mm256_rsqrt_ps(dot256);

  // Newton–Raphson refinement: y = y * (1.5 - 0.5 * dot * y*y)


  __m256 y2 = _mm256_mul_ps(y, y);
  __m256 x_y2 = _mm256_mul_ps(dot256, y2);
  __m256 nr = _mm256_sub_ps(three_halfs, _mm256_mul_ps(half, x_y2));
  y = _mm256_mul_ps(y, nr);

  // Multiply original vector by refined inverse length (normalize)
  __m256 out = _mm256_mul_ps(x, y);

  // Extract x,y,z
  alignas(32) float result[8];
  _mm256_store_ps(result, out);

  return { result[0], result[1], result[2] };
}

#endif

*/


#endif //!VECTOR_UTILITY_H