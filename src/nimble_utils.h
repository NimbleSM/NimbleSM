/*
//@HEADER
// ************************************************************************
//
//                                NimbleSM
//                             Copyright 2018
//   National Technology & Engineering Solutions of Sandia, LLC (NTESS)
//
// Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government
// retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY EXPRESS OR IMPLIED
// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
// MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN
// NO EVENT SHALL NTESS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
// TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions?  Contact David Littlewood (djlittl@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef NIMBLE_UTILS_H
#define NIMBLE_UTILS_H

#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "nimble_defs.h"
#include "nimble_macros.h"
#include "nimble_view.h"

NIMBLE_INLINE_FUNCTION
int
StringLength(const char* str)
{
  int len = 0;
  while (str[len] != '\0') { len++; }
  return len;
}

NIMBLE_INLINE_FUNCTION
bool
StringsAreEqual(const char* str1, const char* str2)
{
  int  len1      = StringLength(str1);
  int  len2      = StringLength(str2);
  int  len       = len1 < len2 ? len1 : len2;
  bool are_equal = true;
  for (int i = 0; i < len; ++i) {
    if (str1[i] != str2[i]) {
      are_equal = false;
      break;
    }
  }
  return are_equal;
}

// Need to avoid conflicts with NimbleSMExtras definitions
#ifndef K_X
static constexpr int K_X = 0;
static constexpr int K_Y = 1;
static constexpr int K_Z = 2;

static constexpr int K_S_XX = 0;
static constexpr int K_S_YY = 1;
static constexpr int K_S_ZZ = 2;
static constexpr int K_S_XY = 3;
static constexpr int K_S_YZ = 4;
static constexpr int K_S_ZX = 5;
static constexpr int K_S_YX = 3;
static constexpr int K_S_ZY = 4;
static constexpr int K_S_XZ = 5;

static constexpr int K_F_XX = 0;
static constexpr int K_F_YY = 1;
static constexpr int K_F_ZZ = 2;
static constexpr int K_F_XY = 3;
static constexpr int K_F_YZ = 4;
static constexpr int K_F_ZX = 5;
static constexpr int K_F_YX = 6;
static constexpr int K_F_ZY = 7;
static constexpr int K_F_XZ = 8;
#endif

template <typename ScalarT>
NIMBLE_INLINE_FUNCTION void
CheckVectorSanity(int vec_length, ScalarT const* const vec, const char* label)
{
#ifdef NIMBLE_DEBUG
  for (int i = 0; i < vec_length; i++) {
#ifndef NIMBLE_HAVE_KOKKOS
    if (!std::isfinite(vec[i])) {
      NIMBLE_ABORT("\n**** Finite value check failed for " + std::string(label) + "!\n");
    }
#else
    if (vec[i] != vec[i] || vec[i] > DBL_MAX || vec[i] < -DBL_MAX) {
      printf("\n**** Error, finite value check failed for %s!\n", label);
    }
#endif
  }
#endif
}

//! Return x times sign of y
template <typename ScalarT>
NIMBLE_INLINE_FUNCTION ScalarT
MultiplySign(const ScalarT x, const ScalarT y)
{
  return x * (1 - 2 * (y < 0));
}

template <typename ScalarT>
NIMBLE_INLINE_FUNCTION ScalarT
Minimum(const ScalarT x, const ScalarT y)
{
  return x < y ? x : y;
}

template <typename ScalarT>
NIMBLE_INLINE_FUNCTION ScalarT
if_then(const bool b, const ScalarT v1, const ScalarT v2)
{
  return b ? v1 : v2;
}

template <typename ScalarT>
NIMBLE_INLINE_FUNCTION ScalarT
if_then_else(const bool b, const ScalarT v1, const ScalarT v2)
{
  return b ? v1 : v2;
}

template <typename ScalarT>
NIMBLE_INLINE_FUNCTION ScalarT
if_then_else_zero(const bool b, const ScalarT v)
{
  return b ? v : ScalarT(0.0);
}

template <typename ScalarT>
void
Print_Full33(std::string label, const ScalarT* const mat, int precision = 2)
{
  std::cout << "Full33 " << label << std::endl;
  std::cout << std::setprecision(precision) << "  " << mat[K_F_XX] << ", " << mat[K_F_XY] << ", " << mat[K_F_XZ]
            << std::endl;
  std::cout << std::setprecision(precision) << "  " << mat[K_F_YX] << ", " << mat[K_F_YY] << ", " << mat[K_F_YZ]
            << std::endl;
  std::cout << std::setprecision(precision) << "  " << mat[K_F_ZX] << ", " << mat[K_F_ZY] << ", " << mat[K_F_ZZ]
            << std::endl;
}

template <typename ScalarT>
void
Print_Sym33(std::string label, const ScalarT* const mat, int precision = 2)
{
  std::cout << "Sym33 " << label << std::endl;
  std::cout << std::setprecision(precision) << "  " << mat[K_S_XX] << ", " << mat[K_S_XY] << ", " << mat[K_S_XZ]
            << std::endl;
  std::cout << std::setprecision(precision) << "  " << mat[K_S_YX] << ", " << mat[K_S_YY] << ", " << mat[K_S_YZ]
            << std::endl;
  std::cout << std::setprecision(precision) << "  " << mat[K_S_ZX] << ", " << mat[K_S_ZY] << ", " << mat[K_S_ZZ]
            << std::endl;
}

//! Compute result = mat^T mat
template <typename ScalarT>
NIMBLE_INLINE_FUNCTION void
Square_Full33T_Full33(const ScalarT* const mat, ScalarT* const result)
{
  NIMBLE_DEBUG_ASSERT(result != mat);
  result[K_S_XX] = mat[K_F_XX] * mat[K_F_XX] + mat[K_F_YX] * mat[K_F_YX] + mat[K_F_ZX] * mat[K_F_ZX];
  result[K_S_YY] = mat[K_F_XY] * mat[K_F_XY] + mat[K_F_YY] * mat[K_F_YY] + mat[K_F_ZY] * mat[K_F_ZY];
  result[K_S_ZZ] = mat[K_F_XZ] * mat[K_F_XZ] + mat[K_F_YZ] * mat[K_F_YZ] + mat[K_F_ZZ] * mat[K_F_ZZ];
  result[K_S_XY] = mat[K_F_XX] * mat[K_F_XY] + mat[K_F_YX] * mat[K_F_YY] + mat[K_F_ZX] * mat[K_F_ZY];
  result[K_S_YZ] = mat[K_F_XY] * mat[K_F_XZ] + mat[K_F_YY] * mat[K_F_YZ] + mat[K_F_ZY] * mat[K_F_ZZ];
  result[K_S_ZX] = mat[K_F_XX] * mat[K_F_XZ] + mat[K_F_YX] * mat[K_F_YZ] + mat[K_F_ZX] * mat[K_F_ZZ];
}

//!  Multiply a full tensor by a scalar
template <typename ScalarT>
NIMBLE_INLINE_FUNCTION void
Mult_Scalar_Full33(ScalarT a, const ScalarT* const mat, ScalarT* const result)
{
  result[K_F_XX] = a * mat[K_F_XX];
  result[K_F_XY] = a * mat[K_F_XY];
  result[K_F_XZ] = a * mat[K_F_XZ];
  result[K_F_YX] = a * mat[K_F_YX];
  result[K_F_YY] = a * mat[K_F_YY];
  result[K_F_YZ] = a * mat[K_F_YZ];
  result[K_F_ZX] = a * mat[K_F_ZX];
  result[K_F_ZY] = a * mat[K_F_ZY];
  result[K_F_ZZ] = a * mat[K_F_ZZ];
}

//!  Sum of a full tensor and a full tensor
template <typename ScalarT>
NIMBLE_INLINE_FUNCTION void
Sum_Full33_Full33(const ScalarT* const A, const ScalarT* const B, ScalarT* const result)
{
  result[K_F_XX] = A[K_F_XX] + B[K_F_XX];
  result[K_F_XY] = A[K_F_XY] + B[K_F_XY];
  result[K_F_XZ] = A[K_F_XZ] + B[K_F_XZ];
  result[K_F_YX] = A[K_F_YX] + B[K_F_YX];
  result[K_F_YY] = A[K_F_YY] + B[K_F_YY];
  result[K_F_YZ] = A[K_F_YZ] + B[K_F_YZ];
  result[K_F_ZX] = A[K_F_ZX] + B[K_F_ZX];
  result[K_F_ZY] = A[K_F_ZY] + B[K_F_ZY];
  result[K_F_ZZ] = A[K_F_ZZ] + B[K_F_ZZ];
}

//!  Sum of a symmetric tensor and a full tensor
template <typename ScalarT>
NIMBLE_INLINE_FUNCTION void
Sum_Sym33_Full33(const ScalarT* const A, const ScalarT* const B, ScalarT* const result)
{
  result[K_F_XX] = A[K_S_XX] + B[K_F_XX];
  result[K_F_XY] = A[K_S_XY] + B[K_F_XY];
  result[K_F_XZ] = A[K_S_XZ] + B[K_F_XZ];
  result[K_F_YX] = A[K_S_YX] + B[K_F_YX];
  result[K_F_YY] = A[K_S_YY] + B[K_F_YY];
  result[K_F_YZ] = A[K_S_YZ] + B[K_F_YZ];
  result[K_F_ZX] = A[K_S_ZX] + B[K_F_ZX];
  result[K_F_ZY] = A[K_S_ZY] + B[K_F_ZY];
  result[K_F_ZZ] = A[K_S_ZZ] + B[K_F_ZZ];
}

//!  Multiply a full tensor by a full tensor
template <typename ScalarT>
NIMBLE_INLINE_FUNCTION void
Mult_Full33_Full33(const ScalarT* const A, const ScalarT* const B, ScalarT* const result)
{
  NIMBLE_DEBUG_ASSERT(A != result);
  NIMBLE_DEBUG_ASSERT(B != result);
  result[K_F_XX] = A[K_F_XX] * B[K_F_XX] + A[K_F_XY] * B[K_F_YX] + A[K_F_XZ] * B[K_F_ZX];
  result[K_F_XY] = A[K_F_XX] * B[K_F_XY] + A[K_F_XY] * B[K_F_YY] + A[K_F_XZ] * B[K_F_ZY];
  result[K_F_XZ] = A[K_F_XX] * B[K_F_XZ] + A[K_F_XY] * B[K_F_YZ] + A[K_F_XZ] * B[K_F_ZZ];
  result[K_F_YX] = A[K_F_YX] * B[K_F_XX] + A[K_F_YY] * B[K_F_YX] + A[K_F_YZ] * B[K_F_ZX];
  result[K_F_YY] = A[K_F_YX] * B[K_F_XY] + A[K_F_YY] * B[K_F_YY] + A[K_F_YZ] * B[K_F_ZY];
  result[K_F_YZ] = A[K_F_YX] * B[K_F_XZ] + A[K_F_YY] * B[K_F_YZ] + A[K_F_YZ] * B[K_F_ZZ];
  result[K_F_ZX] = A[K_F_ZX] * B[K_F_XX] + A[K_F_ZY] * B[K_F_YX] + A[K_F_ZZ] * B[K_F_ZX];
  result[K_F_ZY] = A[K_F_ZX] * B[K_F_XY] + A[K_F_ZY] * B[K_F_YY] + A[K_F_ZZ] * B[K_F_ZY];
  result[K_F_ZZ] = A[K_F_ZX] * B[K_F_XZ] + A[K_F_ZY] * B[K_F_YZ] + A[K_F_ZZ] * B[K_F_ZZ];
}

//!  Multiply a symmetric tensor by a full tensor
template <typename ScalarT>
NIMBLE_INLINE_FUNCTION void
Mult_Sym33_Full33(const ScalarT* const sym, const ScalarT* const full, ScalarT* const result)
{
  NIMBLE_DEBUG_ASSERT(sym != result);
  NIMBLE_DEBUG_ASSERT(full != result);
  result[K_F_XX] = sym[K_S_XX] * full[K_F_XX] + sym[K_S_XY] * full[K_F_YX] + sym[K_S_XZ] * full[K_F_ZX];
  result[K_F_XY] = sym[K_S_XX] * full[K_F_XY] + sym[K_S_XY] * full[K_F_YY] + sym[K_S_XZ] * full[K_F_ZY];
  result[K_F_XZ] = sym[K_S_XX] * full[K_F_XZ] + sym[K_S_XY] * full[K_F_YZ] + sym[K_S_XZ] * full[K_F_ZZ];
  result[K_F_YX] = sym[K_S_YX] * full[K_F_XX] + sym[K_S_YY] * full[K_F_YX] + sym[K_S_YZ] * full[K_F_ZX];
  result[K_F_YY] = sym[K_S_YX] * full[K_F_XY] + sym[K_S_YY] * full[K_F_YY] + sym[K_S_YZ] * full[K_F_ZY];
  result[K_F_YZ] = sym[K_S_YX] * full[K_F_XZ] + sym[K_S_YY] * full[K_F_YZ] + sym[K_S_YZ] * full[K_F_ZZ];
  result[K_F_ZX] = sym[K_S_ZX] * full[K_F_XX] + sym[K_S_ZY] * full[K_F_YX] + sym[K_S_ZZ] * full[K_F_ZX];
  result[K_F_ZY] = sym[K_S_ZX] * full[K_F_XY] + sym[K_S_ZY] * full[K_F_YY] + sym[K_S_ZZ] * full[K_F_ZY];
  result[K_F_ZZ] = sym[K_S_ZX] * full[K_F_XZ] + sym[K_S_ZY] * full[K_F_YZ] + sym[K_S_ZZ] * full[K_F_ZZ];
}

//!  Multiply a full tensor by a full tensor and a scalar
template <typename ScalarT>
NIMBLE_INLINE_FUNCTION void
Mult_Scalar_Full33_Full33(ScalarT alpha, const ScalarT* const A, const ScalarT* const B, ScalarT* const result)
{
  NIMBLE_DEBUG_ASSERT(A != result);
  NIMBLE_DEBUG_ASSERT(B != result);
  result[K_F_XX] = alpha * (A[K_F_XX] * B[K_F_XX] + A[K_F_XY] * B[K_F_YX] + A[K_F_XZ] * B[K_F_ZX]);
  result[K_F_XY] = alpha * (A[K_F_XX] * B[K_F_XY] + A[K_F_XY] * B[K_F_YY] + A[K_F_XZ] * B[K_F_ZY]);
  result[K_F_XZ] = alpha * (A[K_F_XX] * B[K_F_XZ] + A[K_F_XY] * B[K_F_YZ] + A[K_F_XZ] * B[K_F_ZZ]);
  result[K_F_YX] = alpha * (A[K_F_YX] * B[K_F_XX] + A[K_F_YY] * B[K_F_YX] + A[K_F_YZ] * B[K_F_ZX]);
  result[K_F_YY] = alpha * (A[K_F_YX] * B[K_F_XY] + A[K_F_YY] * B[K_F_YY] + A[K_F_YZ] * B[K_F_ZY]);
  result[K_F_YZ] = alpha * (A[K_F_YX] * B[K_F_XZ] + A[K_F_YY] * B[K_F_YZ] + A[K_F_YZ] * B[K_F_ZZ]);
  result[K_F_ZX] = alpha * (A[K_F_ZX] * B[K_F_XX] + A[K_F_ZY] * B[K_F_YX] + A[K_F_ZZ] * B[K_F_ZX]);
  result[K_F_ZY] = alpha * (A[K_F_ZX] * B[K_F_XY] + A[K_F_ZY] * B[K_F_YY] + A[K_F_ZZ] * B[K_F_ZY]);
  result[K_F_ZZ] = alpha * (A[K_F_ZX] * B[K_F_XZ] + A[K_F_ZY] * B[K_F_YZ] + A[K_F_ZZ] * B[K_F_ZZ]);
}

//! Compute result = mat mat^T
template <typename ScalarT>
NIMBLE_INLINE_FUNCTION void
Square_Full33_Full33T(const ScalarT* const mat, ScalarT* const result)
{
  NIMBLE_DEBUG_ASSERT(result != mat);
  result[K_S_XX] = mat[K_F_XX] * mat[K_F_XX] + mat[K_F_XY] * mat[K_F_XY] + mat[K_F_XZ] * mat[K_F_XZ];
  result[K_S_YY] = mat[K_F_YX] * mat[K_F_YX] + mat[K_F_YY] * mat[K_F_YY] + mat[K_F_YZ] * mat[K_F_YZ];
  result[K_S_ZZ] = mat[K_F_ZX] * mat[K_F_ZX] + mat[K_F_ZY] * mat[K_F_ZY] + mat[K_F_ZZ] * mat[K_F_ZZ];
  result[K_S_XY] = mat[K_F_XX] * mat[K_F_YX] + mat[K_F_XY] * mat[K_F_YY] + mat[K_F_XZ] * mat[K_F_YZ];
  result[K_S_YZ] = mat[K_F_YX] * mat[K_F_ZX] + mat[K_F_YY] * mat[K_F_ZY] + mat[K_F_YZ] * mat[K_F_ZZ];
  result[K_S_ZX] = mat[K_F_ZX] * mat[K_F_XX] + mat[K_F_ZY] * mat[K_F_XY] + mat[K_F_ZZ] * mat[K_F_XZ];
}

//!  Multiply a full tensor by a symmetric tensor, return transpose of resultant
template <typename ScalarT>
NIMBLE_INLINE_FUNCTION void
Mult_Full33_Sym33_ReturnT(const ScalarT* const full, const ScalarT* const sym, ScalarT* const result)
{
  NIMBLE_DEBUG_ASSERT(full != result);
  NIMBLE_DEBUG_ASSERT(sym != result);
  result[K_F_XX] = full[K_F_XX] * sym[K_S_XX] + full[K_F_XY] * sym[K_S_YX] + full[K_F_XZ] * sym[K_S_ZX];
  result[K_F_YX] = full[K_F_XX] * sym[K_S_XY] + full[K_F_XY] * sym[K_S_YY] + full[K_F_XZ] * sym[K_S_ZY];
  result[K_F_ZX] = full[K_F_XX] * sym[K_S_XZ] + full[K_F_XY] * sym[K_S_YZ] + full[K_F_XZ] * sym[K_S_ZZ];
  result[K_F_XY] = full[K_F_YX] * sym[K_S_XX] + full[K_F_YY] * sym[K_S_YX] + full[K_F_YZ] * sym[K_S_ZX];
  result[K_F_YY] = full[K_F_YX] * sym[K_S_XY] + full[K_F_YY] * sym[K_S_YY] + full[K_F_YZ] * sym[K_S_ZY];
  result[K_F_ZY] = full[K_F_YX] * sym[K_S_XZ] + full[K_F_YY] * sym[K_S_YZ] + full[K_F_YZ] * sym[K_S_ZZ];
  result[K_F_XZ] = full[K_F_ZX] * sym[K_S_XX] + full[K_F_ZY] * sym[K_S_YX] + full[K_F_ZZ] * sym[K_S_ZX];
  result[K_F_YZ] = full[K_F_ZX] * sym[K_S_XY] + full[K_F_ZY] * sym[K_S_YY] + full[K_F_ZZ] * sym[K_S_ZY];
  result[K_F_ZZ] = full[K_F_ZX] * sym[K_S_XZ] + full[K_F_ZY] * sym[K_S_YZ] + full[K_F_ZZ] * sym[K_S_ZZ];
}

//!  Rotate a symetric tensor, result = R^T S R
template <typename ScalarT>
NIMBLE_INLINE_FUNCTION void
Rotate_Sym33_Using_Rtranspose_S_R(const ScalarT* const s, const ScalarT* const r, ScalarT* const result)
{
  NIMBLE_DEBUG_ASSERT(s != result);
  NIMBLE_DEBUG_ASSERT(r != result);

  ScalarT temp_xx = s[K_S_XX] * r[K_F_XX] + s[K_S_XY] * r[K_F_YX] + s[K_S_XZ] * r[K_F_ZX];
  ScalarT temp_yx = s[K_S_YX] * r[K_F_XX] + s[K_S_YY] * r[K_F_YX] + s[K_S_YZ] * r[K_F_ZX];
  ScalarT temp_zx = s[K_S_ZX] * r[K_F_XX] + s[K_S_ZY] * r[K_F_YX] + s[K_S_ZZ] * r[K_F_ZX];
  ScalarT temp_xy = s[K_S_XX] * r[K_F_XY] + s[K_S_XY] * r[K_F_YY] + s[K_S_XZ] * r[K_F_ZY];
  ScalarT temp_yy = s[K_S_YX] * r[K_F_XY] + s[K_S_YY] * r[K_F_YY] + s[K_S_YZ] * r[K_F_ZY];
  ScalarT temp_zy = s[K_S_ZX] * r[K_F_XY] + s[K_S_ZY] * r[K_F_YY] + s[K_S_ZZ] * r[K_F_ZY];
  ScalarT temp_xz = s[K_S_XX] * r[K_F_XZ] + s[K_S_XY] * r[K_F_YZ] + s[K_S_XZ] * r[K_F_ZZ];
  ScalarT temp_yz = s[K_S_YX] * r[K_F_XZ] + s[K_S_YY] * r[K_F_YZ] + s[K_S_YZ] * r[K_F_ZZ];
  ScalarT temp_zz = s[K_S_ZX] * r[K_F_XZ] + s[K_S_ZY] * r[K_F_YZ] + s[K_S_ZZ] * r[K_F_ZZ];

  result[K_S_XX] = r[K_F_XX] * temp_xx + r[K_F_YX] * temp_yx + r[K_F_ZX] * temp_zx;
  result[K_S_YY] = r[K_F_XY] * temp_xy + r[K_F_YY] * temp_yy + r[K_F_ZY] * temp_zy;
  result[K_S_ZZ] = r[K_F_XZ] * temp_xz + r[K_F_YZ] * temp_yz + r[K_F_ZZ] * temp_zz;
  result[K_S_XY] = r[K_F_XX] * temp_xy + r[K_F_YX] * temp_yy + r[K_F_ZX] * temp_zy;
  result[K_S_YZ] = r[K_F_XY] * temp_xz + r[K_F_YY] * temp_yz + r[K_F_ZY] * temp_zz;
  result[K_S_ZX] = r[K_F_XZ] * temp_xx + r[K_F_YZ] * temp_yx + r[K_F_ZZ] * temp_zx;
}

//!  Unrotate a symetric tensor, result = R S R^T
template <typename ScalarT>
NIMBLE_INLINE_FUNCTION void
Unrotate_Sym33_Using_R_S_Rtranspose(const ScalarT* const s, const ScalarT* const r, ScalarT* const result)
{
  NIMBLE_DEBUG_ASSERT(s != result);
  NIMBLE_DEBUG_ASSERT(r != result);

  ScalarT temp_xx = s[K_S_XX] * r[K_F_XX] + s[K_S_XY] * r[K_F_XY] + s[K_S_XZ] * r[K_F_XZ];
  ScalarT temp_yx = s[K_S_YX] * r[K_F_XX] + s[K_S_YY] * r[K_F_XY] + s[K_S_YZ] * r[K_F_XZ];
  ScalarT temp_zx = s[K_S_ZX] * r[K_F_XX] + s[K_S_ZY] * r[K_F_XY] + s[K_S_ZZ] * r[K_F_XZ];
  ScalarT temp_xy = s[K_S_XX] * r[K_F_YX] + s[K_S_XY] * r[K_F_YY] + s[K_S_XZ] * r[K_F_YZ];
  ScalarT temp_yy = s[K_S_YX] * r[K_F_YX] + s[K_S_YY] * r[K_F_YY] + s[K_S_YZ] * r[K_F_YZ];
  ScalarT temp_zy = s[K_S_ZX] * r[K_F_YX] + s[K_S_ZY] * r[K_F_YY] + s[K_S_ZZ] * r[K_F_YZ];
  ScalarT temp_xz = s[K_S_XX] * r[K_F_ZX] + s[K_S_XY] * r[K_F_ZY] + s[K_S_XZ] * r[K_F_ZZ];
  ScalarT temp_yz = s[K_S_YX] * r[K_F_ZX] + s[K_S_YY] * r[K_F_ZY] + s[K_S_YZ] * r[K_F_ZZ];
  ScalarT temp_zz = s[K_S_ZX] * r[K_F_ZX] + s[K_S_ZY] * r[K_F_ZY] + s[K_S_ZZ] * r[K_F_ZZ];

  result[K_S_XX] = r[K_F_XX] * temp_xx + r[K_F_XY] * temp_yx + r[K_F_XZ] * temp_zx;
  result[K_S_YY] = r[K_F_YX] * temp_xy + r[K_F_YY] * temp_yy + r[K_F_YZ] * temp_zy;
  result[K_S_ZZ] = r[K_F_ZX] * temp_xz + r[K_F_ZY] * temp_yz + r[K_F_ZZ] * temp_zz;
  result[K_S_XY] = r[K_F_XX] * temp_xy + r[K_F_XY] * temp_yy + r[K_F_XZ] * temp_zy;
  result[K_S_YZ] = r[K_F_YX] * temp_xz + r[K_F_YY] * temp_yz + r[K_F_YZ] * temp_zz;
  result[K_S_ZX] = r[K_F_ZX] * temp_xx + r[K_F_ZY] * temp_yx + r[K_F_ZZ] * temp_zx;
}

template <typename ScalarT>
NIMBLE_INLINE_FUNCTION void
CrossProduct(const ScalarT* const u, const ScalarT* const v, ScalarT* const result)
{
  result[0] = u[1] * v[2] - u[2] * v[1];
  result[1] = u[2] * v[0] - u[0] * v[2];
  result[2] = u[0] * v[1] - u[1] * v[0];
}

template <typename ScalarT>
NIMBLE_INLINE_FUNCTION ScalarT
Determinant_Full33(const ScalarT* const mat)
{
  ScalarT minor0 = mat[K_F_YY] * mat[K_F_ZZ] - mat[K_F_YZ] * mat[K_F_ZY];
  ScalarT minor1 = mat[K_F_YX] * mat[K_F_ZZ] - mat[K_F_YZ] * mat[K_F_ZX];
  ScalarT minor2 = mat[K_F_YX] * mat[K_F_ZY] - mat[K_F_YY] * mat[K_F_ZX];
  ScalarT det    = mat[K_F_XX] * minor0 - mat[K_F_XY] * minor1 + mat[K_F_XZ] * minor2;

  return det;
}

template <typename ScalarT>
NIMBLE_INLINE_FUNCTION ScalarT
Norm_Sym33(const ScalarT* const mat)
{
  ScalarT norm = 0.0;
  norm += mat[K_S_XX] * mat[K_S_XX] + mat[K_S_XY] * mat[K_S_XY] + mat[K_S_XZ] * mat[K_S_XZ];
  norm += mat[K_S_YX] * mat[K_S_YX] + mat[K_S_YY] * mat[K_S_YY] + mat[K_S_YZ] * mat[K_S_YZ];
  norm += mat[K_S_ZX] * mat[K_S_ZX] + mat[K_S_ZY] * mat[K_S_ZY] + mat[K_S_ZZ] * mat[K_S_ZZ];

  if (norm > 0.0) { norm = std::sqrt(norm); }

  return norm;
}

template <typename ScalarT>
NIMBLE_INLINE_FUNCTION ScalarT
Norm_Full33(const ScalarT* const mat)
{
  ScalarT norm = 0.0;
  norm += mat[K_F_XX] * mat[K_F_XX] + mat[K_F_XY] * mat[K_F_XY] + mat[K_F_XZ] * mat[K_F_XZ];
  norm += mat[K_F_YX] * mat[K_F_YX] + mat[K_F_YY] * mat[K_F_YY] + mat[K_F_YZ] * mat[K_F_YZ];
  norm += mat[K_F_ZX] * mat[K_F_ZX] + mat[K_F_ZY] * mat[K_F_ZY] + mat[K_F_ZZ] * mat[K_F_ZZ];

  if (norm > 0.0) { norm = std::sqrt(norm); }

  return norm;
}

template <typename ScalarT>
NIMBLE_INLINE_FUNCTION void
Zero_Full33(ScalarT* const mat)
{
  mat[K_F_XX] = 0.0;
  mat[K_F_XY] = 0.0;
  mat[K_F_XZ] = 0.0;
  mat[K_F_YX] = 0.0;
  mat[K_F_YY] = 0.0;
  mat[K_F_YZ] = 0.0;
  mat[K_F_ZX] = 0.0;
  mat[K_F_ZY] = 0.0;
  mat[K_F_ZZ] = 0.0;
}

template <typename ScalarT>
NIMBLE_INLINE_FUNCTION void
SetEqual_Full33(const ScalarT* const B, ScalarT* const A)
{
  A[K_F_XX] = B[K_F_XX];
  A[K_F_XY] = B[K_F_XY];
  A[K_F_XZ] = B[K_F_XZ];
  A[K_F_YX] = B[K_F_YX];
  A[K_F_YY] = B[K_F_YY];
  A[K_F_YZ] = B[K_F_YZ];
  A[K_F_ZX] = B[K_F_ZX];
  A[K_F_ZY] = B[K_F_ZY];
  A[K_F_ZZ] = B[K_F_ZZ];
}

template <typename ScalarT>
NIMBLE_INLINE_FUNCTION void
SetEqual_Sym33(const ScalarT* const B, ScalarT* const A)
{
  A[K_S_XX] = B[K_S_XX];
  A[K_S_YY] = B[K_S_YY];
  A[K_S_ZZ] = B[K_S_ZZ];
  A[K_S_XY] = B[K_S_XY];
  A[K_S_YZ] = B[K_S_YZ];
  A[K_S_ZX] = B[K_S_ZX];
}

template <typename ScalarT>
NIMBLE_INLINE_FUNCTION void
SymPart_Full33(const ScalarT* const mat, ScalarT* const result)
{
  result[K_S_XX] = mat[K_F_XX];
  result[K_S_YY] = mat[K_F_YY];
  result[K_S_ZZ] = mat[K_F_ZZ];
  result[K_S_XY] = 0.5 * (mat[K_F_XY] + mat[K_F_YX]);
  result[K_S_YZ] = 0.5 * (mat[K_F_YZ] + mat[K_F_ZY]);
  result[K_S_ZX] = 0.5 * (mat[K_F_ZX] + mat[K_F_XZ]);
}

template <typename ScalarT>
NIMBLE_INLINE_FUNCTION void
SkewPart_Full33(const ScalarT* const mat, ScalarT* const result)
{
  result[K_F_XX] = 0.0;
  result[K_F_XY] = 0.5 * (mat[K_F_XY] - mat[K_F_YX]);
  result[K_F_XZ] = 0.5 * (mat[K_F_XZ] - mat[K_F_ZX]);
  result[K_F_YX] = 0.5 * (mat[K_F_YX] - mat[K_F_XY]);
  result[K_F_YY] = 0.0;
  result[K_F_YZ] = 0.5 * (mat[K_F_YZ] - mat[K_F_ZY]);
  result[K_F_ZX] = 0.5 * (mat[K_F_ZX] - mat[K_F_XZ]);
  result[K_F_ZY] = 0.5 * (mat[K_F_ZY] - mat[K_F_YZ]);
  result[K_F_ZZ] = 0.0;
}

/// \brief Function to invert a 3x3 matrix
///
/// \tparam ScalarT  Scalar values
/// \param mat  Pointer to 3x3 matrix with a storage based on the indices K_F_**
/// \param inv  Pointer to 3x3 mqtrix with a storage based on the indices K_F_**
/// \return  Determinant
/// \remark  The matrix is has the pattern
///           [  K_F_XX  K_F_XY  K_F_XZ  ]
///           [  K_F_YX  K_F_YY  K_F_YZ  ]
///           [  K_F_ZX  K_F_ZY  K_F_ZZ  ]
/// but note that K_F_YY = 1 and K_F_ZZ = 2
/// The storage proceeds along diagonals.
template <typename ScalarT>
NIMBLE_INLINE_FUNCTION ScalarT
Invert_Full33(const ScalarT* mat, ScalarT* inv)
{
  ScalarT minor0 = mat[K_F_YY] * mat[K_F_ZZ] - mat[K_F_YZ] * mat[K_F_ZY];
  ScalarT minor1 = mat[K_F_YX] * mat[K_F_ZZ] - mat[K_F_YZ] * mat[K_F_ZX];
  ScalarT minor2 = mat[K_F_YX] * mat[K_F_ZY] - mat[K_F_YY] * mat[K_F_ZX];
  ScalarT minor3 = mat[K_F_XY] * mat[K_F_ZZ] - mat[K_F_XZ] * mat[K_F_ZY];
  ScalarT minor4 = mat[K_F_XX] * mat[K_F_ZZ] - mat[K_F_ZX] * mat[K_F_XZ];
  ScalarT minor5 = mat[K_F_XX] * mat[K_F_ZY] - mat[K_F_XY] * mat[K_F_ZX];
  ScalarT minor6 = mat[K_F_XY] * mat[K_F_YZ] - mat[K_F_XZ] * mat[K_F_YY];
  ScalarT minor7 = mat[K_F_XX] * mat[K_F_YZ] - mat[K_F_XZ] * mat[K_F_YX];
  ScalarT minor8 = mat[K_F_XX] * mat[K_F_YY] - mat[K_F_XY] * mat[K_F_YX];
  ScalarT det    = mat[K_F_XX] * minor0 - mat[K_F_XY] * minor1 + mat[K_F_XZ] * minor2;

  NIMBLE_ASSERT(det >= 0.0, "Error in Invert_Full33(), singular matrix");

  inv[K_F_XX] = minor0 / det;
  inv[K_F_XY] = -1.0 * minor3 / det;
  inv[K_F_XZ] = minor6 / det;
  inv[K_F_YX] = -1.0 * minor1 / det;
  inv[K_F_YY] = minor4 / det;
  inv[K_F_YZ] = -1.0 * minor7 / det;
  inv[K_F_ZX] = minor2 / det;
  inv[K_F_ZY] = -1.0 * minor5 / det;
  inv[K_F_ZZ] = minor8 / det;

  return det;
}

//! R^N logarithmic map using Baker-Campbell-Hausdorff (BCH) expansion (4 terms)
// All arguments are assumed to be full tensors (9 components)
template <typename ScalarT>
NIMBLE_INLINE_FUNCTION void
BCH(const ScalarT* const x, const ScalarT* const y, ScalarT* const result)
{
  // adapted from MiniTensor_LinearAlgebra.t.h bch()

  ScalarT temp1[9], temp2[9], temp3[9];

  // first order term
  // x + y
  Sum_Full33_Full33(x, y, result);

  // second order term
  // 0.5*(x*y - y*x)
  Mult_Scalar_Full33_Full33(0.5, x, y, temp1);
  Mult_Scalar_Full33_Full33(-0.5, y, x, temp2);
  Sum_Full33_Full33(result, temp1, result);
  Sum_Full33_Full33(result, temp2, result);

  // third order term
  // 1.0/12.0 * (x*x*y - 2.0*x*y*x + x*y*y + y*x*x - 2.0*y*x*y + y*y*x)
  // 1/12 * x*x*y
  Mult_Full33_Full33(x, y, temp1);
  Mult_Scalar_Full33_Full33(1.0 / 12.0, x, temp1, temp2);
  Sum_Full33_Full33(result, temp2, result);
  // -1/6 * x*y*x
  Mult_Full33_Full33(y, x, temp1);
  Mult_Scalar_Full33_Full33(-1.0 / 6.0, x, temp1, temp2);
  Sum_Full33_Full33(result, temp2, result);
  // 1/12 * x*y*y
  Mult_Full33_Full33(y, y, temp1);
  Mult_Scalar_Full33_Full33(1.0 / 12.0, x, temp1, temp2);
  Sum_Full33_Full33(result, temp2, result);
  // 1/12 * y*x*x
  Mult_Full33_Full33(x, x, temp1);
  Mult_Scalar_Full33_Full33(1.0 / 12.0, y, temp1, temp2);
  Sum_Full33_Full33(result, temp2, result);
  // -1/6 * y*x*y
  Mult_Full33_Full33(x, y, temp1);
  Mult_Scalar_Full33_Full33(-1.0 / 6.0, y, temp1, temp2);
  Sum_Full33_Full33(result, temp2, result);
  // 1/12 * y*y*x
  Mult_Full33_Full33(y, x, temp1);
  Mult_Scalar_Full33_Full33(1.0 / 12.0, y, temp1, temp2);
  Sum_Full33_Full33(result, temp2, result);

  // fourth order term
  // 1.0/24.0 * (x*x*y*y - 2.0*x*y*x*y + 2.0*y*x*y*x - y*y*x*x);
  // 1.0/24.0 * x*x*y*y
  Mult_Full33_Full33(y, y, temp1);
  Mult_Full33_Full33(x, temp1, temp2);
  Mult_Scalar_Full33_Full33(1.0 / 24.0, x, temp2, temp3);
  Sum_Full33_Full33(result, temp3, result);
  // -2.0/24.0 * x*y*x*y
  Mult_Full33_Full33(x, y, temp1);
  Mult_Full33_Full33(y, temp1, temp2);
  Mult_Scalar_Full33_Full33(-2.0 / 24.0, x, temp2, temp3);
  Sum_Full33_Full33(result, temp3, result);
  // 2.0/24.0 * y*x*y*x
  Mult_Full33_Full33(y, x, temp1);
  Mult_Full33_Full33(x, temp1, temp2);
  Mult_Scalar_Full33_Full33(2.0 / 24.0, y, temp2, temp3);
  Sum_Full33_Full33(result, temp3, result);
  // -1.0/24.0 * y*y*x*x
  Mult_Full33_Full33(x, x, temp1);
  Mult_Full33_Full33(y, temp1, temp2);
  Mult_Scalar_Full33_Full33(-1.0 / 24.0, y, temp2, temp3);
  Sum_Full33_Full33(result, temp3, result);
}

//! R^N logarithmic map using Baker-Campbell-Hausdorff (BCH) expansion (4 terms)
// All arguments are assumed to be full tensors (9 components)
template <typename ScalarT>
NIMBLE_INLINE_FUNCTION void
BCH_Sym33_Full33(const ScalarT* const sym, const ScalarT* const full, ScalarT* const result)
{
  ScalarT temp[9];
  temp[K_F_XX] = sym[K_S_XX];
  temp[K_F_XY] = sym[K_S_XY];
  temp[K_F_XZ] = sym[K_S_XZ];
  temp[K_F_YX] = sym[K_S_YX];
  temp[K_F_YY] = sym[K_S_YY];
  temp[K_F_YZ] = sym[K_S_YZ];
  temp[K_F_ZX] = sym[K_S_ZX];
  temp[K_F_ZY] = sym[K_S_ZY];
  temp[K_F_ZZ] = sym[K_S_ZZ];
  BCH(temp, full, result);
}

//! Pade approximation to cos( acos(x)/3 )
// Mathematica commands:
// Needs["FunctionApproximations`"]
// r1 = MiniMaxApproximation[Cos[ArcCos[x]/3], {x, {0, 1}, 6, 5},
// WorkingPrecision -> 18, MaxIterations -> 500] 6 and 5 indicate the polynomial
// order in the numerator and denominator.
template <typename ScalarT>
NIMBLE_INLINE_FUNCTION ScalarT
Cos_Of_Acos_Divided_By_3(const ScalarT x)
{
  // algorithm adapted from Math_Trig.h

  const ScalarT x2 = x * x;
  const ScalarT x4 = x2 * x2;

  return (0.866025403784438713 + 2.12714890259493060 * x +
          ((1.89202064815951569 + 0.739603278343401613 * x) * x2 +
           (0.121973926953064794 + x * (0.00655637626263929360 + 0.0000390884982780803443 * x)) * x4)) /
         (1.0 + 2.26376989330935617 * x +
          ((1.80461009751278976 + 0.603976798217196003 * x) * x2 +
           (0.0783255761115461708 + 0.00268525944538021629 * x) * x4));
}

template <typename ScalarT>
NIMBLE_INLINE_FUNCTION void
Eigen_Sym33_NonUnit(
    const ScalarT* const tensor,
    ScalarT&             eval0,
    ScalarT&             eval1,
    ScalarT&             eval2,
    ScalarT* const       evec0,
    ScalarT* const       evec1,
    ScalarT* const       evec2)
{
  ScalarT        cxx = tensor[K_S_XX];
  ScalarT        cyy = tensor[K_S_YY];
  ScalarT        czz = tensor[K_S_ZZ];
  const ScalarT& cxy = tensor[K_S_XY];
  const ScalarT& cyz = tensor[K_S_YZ];
  const ScalarT& czx = tensor[K_S_ZX];

  const ScalarT c1 = (cxx + cyy + czz) / ScalarT(3.0);

  cxx -= c1;
  cyy -= c1;
  czz -= c1;

  const ScalarT cxy_cxy = cxy * cxy;
  const ScalarT cyz_cyz = cyz * cyz;
  const ScalarT czx_czx = czx * czx;
  const ScalarT cxx_cyy = cxx * cyy;

  const ScalarT c2 = cxx_cyy + cyy * czz + czz * cxx - cxy_cxy - cyz_cyz - czx_czx;

  const ScalarT ThreeOverA     = ScalarT(-3.0) / c2;
  const ScalarT sqrtThreeOverA = std::sqrt(ThreeOverA);

  const ScalarT c3 = cxx * cyz_cyz + cyy * czx_czx - ScalarT(2.0) * cxy * cyz * czx + czz * (cxy_cxy - cxx_cyy);

  const ScalarT rr = ScalarT(-0.5) * c3 * ThreeOverA * sqrtThreeOverA;

  const ScalarT arg = Minimum(std::abs(rr), ScalarT(1.0));

  const ScalarT cos_thd3 = Cos_Of_Acos_Divided_By_3(arg);

  const ScalarT two_cos_thd3 = ScalarT(2.0) * MultiplySign(cos_thd3, rr);

  eval2 = two_cos_thd3 / sqrtThreeOverA;

  ScalarT crow0[3] = {cxx - eval2, cxy, czx};
  ScalarT crow1[3] = {cxy, cyy - eval2, cyz};
  ScalarT crow2[3] = {czx, cyz, czz - eval2};

  //
  // do QR decomposition with column pivoting
  //
  const ScalarT k0 = crow0[0] * crow0[0] + cxy_cxy + czx_czx;
  const ScalarT k1 = cxy_cxy + crow1[1] * crow1[1] + cyz_cyz;
  const ScalarT k2 = czx_czx + cyz_cyz + crow2[2] * crow2[2];

  ScalarT k_row1[3];
  ScalarT row2[3];
  ScalarT row3[3];

  // returns zero or nan
  const bool k0gk1 = k1 <= k0;
  const bool k0gk2 = k2 <= k0;
  const bool k1gk2 = k2 <= k1;

  const bool k0_largest = k0gk1 && k0gk2;
  const bool k1_largest = k1gk2 && !k0gk1;
  const bool k2_largest = !(k0_largest || k1_largest);

  k_row1[0] = if_then_else_zero(k0_largest, crow0[0]) + if_then_else_zero(k1_largest, crow1[0]) +
              if_then_else_zero(k2_largest, crow2[0]);

  k_row1[1] = if_then_else_zero(k0_largest, crow0[1]) + if_then_else_zero(k1_largest, crow1[1]) +
              if_then_else_zero(k2_largest, crow2[1]);

  k_row1[2] = if_then_else_zero(k0_largest, crow0[2]) + if_then_else_zero(k1_largest, crow1[2]) +
              if_then_else_zero(k2_largest, crow2[2]);

  row2[0] = if_then_else(k0_largest, crow1[0], crow0[0]);
  row2[1] = if_then_else(k0_largest, crow1[1], crow0[1]);
  row2[2] = if_then_else(k0_largest, crow1[2], crow0[2]);

  row3[0] = if_then_else(k2_largest, crow1[0], crow2[0]);
  row3[1] = if_then_else(k2_largest, crow1[1], crow2[1]);
  row3[2] = if_then_else(k2_largest, crow1[2], crow2[2]);

  const ScalarT ki_ki = ScalarT(1.0) / (if_then_else_zero(k0_largest, k0) + if_then_else_zero(k1_largest, k1) +
                                        if_then_else_zero(k2_largest, k2));

  const ScalarT ki_dpr1 = ki_ki * (k_row1[0] * row2[0] + k_row1[1] * row2[1] + k_row1[2] * row2[2]);
  const ScalarT ki_dpr2 = ki_ki * (k_row1[0] * row3[0] + k_row1[1] * row3[1] + k_row1[2] * row3[2]);

  row2[0] -= ki_dpr1 * k_row1[0];
  row2[1] -= ki_dpr1 * k_row1[1];
  row2[2] -= ki_dpr1 * k_row1[2];

  row3[0] -= ki_dpr2 * k_row1[0];
  row3[1] -= ki_dpr2 * k_row1[1];
  row3[2] -= ki_dpr2 * k_row1[2];

  const ScalarT a0 = row2[0] * row2[0] + row2[1] * row2[1] + row2[2] * row2[2];
  const ScalarT a1 = row3[0] * row3[0] + row3[1] * row3[1] + row3[2] * row3[2];

  ScalarT a_row2[3];

  const bool a0lea1 = a0 <= a1;

  a_row2[0]           = if_then_else(a0lea1, row3[0], row2[0]);
  a_row2[1]           = if_then_else(a0lea1, row3[1], row2[1]);
  a_row2[2]           = if_then_else(a0lea1, row3[2], row2[2]);
  const ScalarT ai_ai = ScalarT(1.0) / if_then_else(a0lea1, a1, a0);

  evec2[K_X] = k_row1[1] * a_row2[2] - k_row1[2] * a_row2[1];
  evec2[K_Y] = k_row1[2] * a_row2[0] - k_row1[0] * a_row2[2];
  evec2[K_Z] = k_row1[0] * a_row2[1] - k_row1[1] * a_row2[0];

  const ScalarT k_atr11 = cxx * k_row1[0] + cxy * k_row1[1] + czx * k_row1[2];
  const ScalarT k_atr21 = cxy * k_row1[0] + cyy * k_row1[1] + cyz * k_row1[2];
  const ScalarT k_atr31 = czx * k_row1[0] + cyz * k_row1[1] + czz * k_row1[2];

  const ScalarT a_atr12 = cxx * a_row2[0] + cxy * a_row2[1] + czx * a_row2[2];
  const ScalarT a_atr22 = cxy * a_row2[0] + cyy * a_row2[1] + cyz * a_row2[2];
  const ScalarT a_atr32 = czx * a_row2[0] + cyz * a_row2[1] + czz * a_row2[2];

  ScalarT       rm2xx       = (k_row1[0] * k_atr11 + k_row1[1] * k_atr21 + k_row1[2] * k_atr31) * ki_ki;
  const ScalarT k_a_rm2xy   = (k_row1[0] * a_atr12 + k_row1[1] * a_atr22 + k_row1[2] * a_atr32);
  ScalarT       rm2yy       = (a_row2[0] * a_atr12 + a_row2[1] * a_atr22 + a_row2[2] * a_atr32) * ai_ai;
  const ScalarT rm2xy_rm2xy = k_a_rm2xy * k_a_rm2xy * ai_ai * ki_ki;

  //
  // Wilkinson shift
  //
  const ScalarT b        = ScalarT(0.5) * (rm2xx - rm2yy);
  const ScalarT sqrtTerm = MultiplySign(std::sqrt(b * b + rm2xy_rm2xy), b);
  eval0                  = rm2yy + b - sqrtTerm;

  eval1 = rm2xx + rm2yy - eval0;

  rm2xx -= eval0;
  rm2yy -= eval0;

  const ScalarT rm2xx2 = rm2xx * rm2xx;
  const ScalarT rm2yy2 = rm2yy * rm2yy;

  const ScalarT fac1 = if_then_else(rm2xx2 < rm2yy2, k_a_rm2xy * ai_ai, rm2xx);
  const ScalarT fac2 = if_then_else(rm2xx2 < rm2yy2, rm2yy, ki_ki * k_a_rm2xy);

  evec0[0] = fac1 * a_row2[0] - fac2 * k_row1[0];
  evec0[1] = fac1 * a_row2[1] - fac2 * k_row1[1];
  evec0[2] = fac1 * a_row2[2] - fac2 * k_row1[2];

  const bool rm2xx2iszero      = rm2xx2 == ScalarT(0.0);
  const bool rm2xy_rm2xyiszero = rm2xy_rm2xy == ScalarT(0.0);
  const bool both_zero         = rm2xx2iszero && rm2xy_rm2xyiszero;

  // check degeneracy
  evec0[0] = if_then_else(both_zero, a_row2[0], evec0[0]);
  evec0[1] = if_then_else(both_zero, a_row2[1], evec0[1]);
  evec0[2] = if_then_else(both_zero, a_row2[2], evec0[2]);

  evec1[K_X] = evec2[K_Y] * evec0[K_Z] - evec2[K_Z] * evec0[K_Y];
  evec1[K_Y] = evec2[K_Z] * evec0[K_X] - evec2[K_X] * evec0[K_Z];
  evec1[K_Z] = evec2[K_X] * evec0[K_Y] - evec2[K_Y] * evec0[K_X];

  eval0 += c1;
  eval1 += c1;
  eval2 += c1;

  const ScalarT c2tol = (c1 * c1) * (-1.0e-30);  // DJL where does the magic number come from?

  const bool c2lsmall_neg = c2 < c2tol;

  eval0 = if_then_else(c2lsmall_neg, eval0, c1);
  eval1 = if_then_else(c2lsmall_neg, eval1, c1);
  eval2 = if_then_else(c2lsmall_neg, eval2, c1);

  evec0[0] = if_then_else(c2lsmall_neg, evec0[0], ScalarT(1.0));
  evec0[1] = if_then_else(c2lsmall_neg, evec0[1], ScalarT(0.0));
  evec0[2] = if_then_else(c2lsmall_neg, evec0[2], ScalarT(0.0));

  evec1[0] = if_then_else(c2lsmall_neg, evec1[0], ScalarT(0.0));
  evec1[1] = if_then_else(c2lsmall_neg, evec1[1], ScalarT(1.0));
  evec1[2] = if_then_else(c2lsmall_neg, evec1[2], ScalarT(0.0));

  evec2[0] = if_then_else(c2lsmall_neg, evec2[0], ScalarT(0.0));
  evec2[1] = if_then_else(c2lsmall_neg, evec2[1], ScalarT(0.0));
  evec2[2] = if_then_else(c2lsmall_neg, evec2[2], ScalarT(1.0));

  return;
}

template <typename ScalarT>
NIMBLE_INLINE_FUNCTION void
Eigen_Sym33(
    const ScalarT* const tensor,
    ScalarT&             eval0,
    ScalarT&             eval1,
    ScalarT&             eval2,
    ScalarT* const       evec0,
    ScalarT* const       evec1,
    ScalarT* const       evec2)
{
  ScalarT        cxx = tensor[K_S_XX];
  ScalarT        cyy = tensor[K_S_YY];
  ScalarT        czz = tensor[K_S_ZZ];
  const ScalarT& cxy = tensor[K_S_XY];
  const ScalarT& cyz = tensor[K_S_YZ];
  const ScalarT& czx = tensor[K_S_ZX];

  const ScalarT c1 = (cxx + cyy + czz) / ScalarT(3.0);

  cxx -= c1;
  cyy -= c1;
  czz -= c1;

  const ScalarT cxy_cxy = cxy * cxy;
  const ScalarT cyz_cyz = cyz * cyz;
  const ScalarT czx_czx = czx * czx;
  const ScalarT cxx_cyy = cxx * cyy;

  const ScalarT c2 = cxx_cyy + cyy * czz + czz * cxx - cxy_cxy - cyz_cyz - czx_czx;
  const ScalarT c3 = cxx * cyz_cyz + cyy * czx_czx - ScalarT(2.0) * cxy * cyz * czx + czz * (cxy_cxy - cxx_cyy);

  // initialize eigenvalues and eigenvectors to those for the identity matrix
  eval0 = c1;
  eval1 = c1;
  eval2 = c1;

  evec0[0] = ScalarT(1.0);
  evec0[1] = ScalarT(0.0);
  evec0[2] = ScalarT(0.0);

  evec1[0] = ScalarT(0.0);
  evec1[1] = ScalarT(1.0);
  evec1[2] = ScalarT(0.0);

  evec2[0] = ScalarT(0.0);
  evec2[1] = ScalarT(0.0);
  evec2[2] = ScalarT(1.0);

  // if c2 is near zero, then tensor is identity (?); no further comutation needed
  const ScalarT c2tol = (c1 * c1) * (-1.0e-30);  // DJL where does the magic number come from?
  const bool c2lsmall_neg = c2 < c2tol;
  if (!c2lsmall_neg) {
      return;
  }

  const ScalarT ThreeOverA     = ScalarT(-3.0) / c2;
  const ScalarT sqrtThreeOverA = std::sqrt(ThreeOverA);

  const ScalarT rr = ScalarT(-0.5) * c3 * ThreeOverA * sqrtThreeOverA;

  const ScalarT arg = Minimum(std::abs(rr), ScalarT(1.0));

  const ScalarT cos_thd3 = Cos_Of_Acos_Divided_By_3(arg);

  const ScalarT two_cos_thd3 = ScalarT(2.0) * MultiplySign(cos_thd3, rr);

  eval2 = two_cos_thd3 / sqrtThreeOverA;

  ScalarT crow0[3] = {cxx - eval2, cxy, czx};
  ScalarT crow1[3] = {cxy, cyy - eval2, cyz};
  ScalarT crow2[3] = {czx, cyz, czz - eval2};

  //
  // do QR decomposition with column pivoting
  //
  const ScalarT k0 = crow0[0] * crow0[0] + cxy_cxy + czx_czx;
  const ScalarT k1 = cxy_cxy + crow1[1] * crow1[1] + cyz_cyz;
  const ScalarT k2 = czx_czx + cyz_cyz + crow2[2] * crow2[2];

  ScalarT k_row1[3];
  ScalarT row2[3];
  ScalarT row3[3];

  // returns zero or nan
  const bool k0gk1 = k1 <= k0;
  const bool k0gk2 = k2 <= k0;
  const bool k1gk2 = k2 <= k1;

  const bool k0_largest = k0gk1 && k0gk2;
  const bool k1_largest = k1gk2 && !k0gk1;
  const bool k2_largest = !(k0_largest || k1_largest);

  k_row1[0] = if_then_else_zero(k0_largest, crow0[0]) + if_then_else_zero(k1_largest, crow1[0]) +
              if_then_else_zero(k2_largest, crow2[0]);

  k_row1[1] = if_then_else_zero(k0_largest, crow0[1]) + if_then_else_zero(k1_largest, crow1[1]) +
              if_then_else_zero(k2_largest, crow2[1]);

  k_row1[2] = if_then_else_zero(k0_largest, crow0[2]) + if_then_else_zero(k1_largest, crow1[2]) +
              if_then_else_zero(k2_largest, crow2[2]);

  row2[0] = if_then_else(k0_largest, crow1[0], crow0[0]);
  row2[1] = if_then_else(k0_largest, crow1[1], crow0[1]);
  row2[2] = if_then_else(k0_largest, crow1[2], crow0[2]);

  row3[0] = if_then_else(k2_largest, crow1[0], crow2[0]);
  row3[1] = if_then_else(k2_largest, crow1[1], crow2[1]);
  row3[2] = if_then_else(k2_largest, crow1[2], crow2[2]);

  const ScalarT ki_ki = ScalarT(1.0) / (if_then_else_zero(k0_largest, k0) + if_then_else_zero(k1_largest, k1) +
                                        if_then_else_zero(k2_largest, k2));

  const ScalarT ki_dpr1 = ki_ki * (k_row1[0] * row2[0] + k_row1[1] * row2[1] + k_row1[2] * row2[2]);
  const ScalarT ki_dpr2 = ki_ki * (k_row1[0] * row3[0] + k_row1[1] * row3[1] + k_row1[2] * row3[2]);

  row2[0] -= ki_dpr1 * k_row1[0];
  row2[1] -= ki_dpr1 * k_row1[1];
  row2[2] -= ki_dpr1 * k_row1[2];

  row3[0] -= ki_dpr2 * k_row1[0];
  row3[1] -= ki_dpr2 * k_row1[1];
  row3[2] -= ki_dpr2 * k_row1[2];

  const ScalarT a0 = row2[0] * row2[0] + row2[1] * row2[1] + row2[2] * row2[2];
  const ScalarT a1 = row3[0] * row3[0] + row3[1] * row3[1] + row3[2] * row3[2];

  ScalarT a_row2[3];

  const bool a0lea1 = a0 <= a1;

  a_row2[0]           = if_then_else(a0lea1, row3[0], row2[0]);
  a_row2[1]           = if_then_else(a0lea1, row3[1], row2[1]);
  a_row2[2]           = if_then_else(a0lea1, row3[2], row2[2]);
  const ScalarT ai_ai = ScalarT(1.0) / if_then_else(a0lea1, a1, a0);

  evec2[K_X] = k_row1[1] * a_row2[2] - k_row1[2] * a_row2[1];
  evec2[K_Y] = k_row1[2] * a_row2[0] - k_row1[0] * a_row2[2];
  evec2[K_Z] = k_row1[0] * a_row2[1] - k_row1[1] * a_row2[0];

  const ScalarT k_atr11 = cxx * k_row1[0] + cxy * k_row1[1] + czx * k_row1[2];
  const ScalarT k_atr21 = cxy * k_row1[0] + cyy * k_row1[1] + cyz * k_row1[2];
  const ScalarT k_atr31 = czx * k_row1[0] + cyz * k_row1[1] + czz * k_row1[2];

  const ScalarT a_atr12 = cxx * a_row2[0] + cxy * a_row2[1] + czx * a_row2[2];
  const ScalarT a_atr22 = cxy * a_row2[0] + cyy * a_row2[1] + cyz * a_row2[2];
  const ScalarT a_atr32 = czx * a_row2[0] + cyz * a_row2[1] + czz * a_row2[2];

  ScalarT       rm2xx       = (k_row1[0] * k_atr11 + k_row1[1] * k_atr21 + k_row1[2] * k_atr31) * ki_ki;
  const ScalarT k_a_rm2xy   = (k_row1[0] * a_atr12 + k_row1[1] * a_atr22 + k_row1[2] * a_atr32);
  ScalarT       rm2yy       = (a_row2[0] * a_atr12 + a_row2[1] * a_atr22 + a_row2[2] * a_atr32) * ai_ai;
  const ScalarT rm2xy_rm2xy = k_a_rm2xy * k_a_rm2xy * ai_ai * ki_ki;

  //
  // Wilkinson shift
  //
  const ScalarT b        = ScalarT(0.5) * (rm2xx - rm2yy);
  const ScalarT sqrtTerm = MultiplySign(std::sqrt(b * b + rm2xy_rm2xy), b);
  eval0                  = rm2yy + b - sqrtTerm;

  eval1 = rm2xx + rm2yy - eval0;

  rm2xx -= eval0;
  rm2yy -= eval0;

  const ScalarT rm2xx2 = rm2xx * rm2xx;
  const ScalarT rm2yy2 = rm2yy * rm2yy;

  const ScalarT fac1 = if_then_else(rm2xx2 < rm2yy2, k_a_rm2xy * ai_ai, rm2xx);
  const ScalarT fac2 = if_then_else(rm2xx2 < rm2yy2, rm2yy, ki_ki * k_a_rm2xy);

  evec0[0] = fac1 * a_row2[0] - fac2 * k_row1[0];
  evec0[1] = fac1 * a_row2[1] - fac2 * k_row1[1];
  evec0[2] = fac1 * a_row2[2] - fac2 * k_row1[2];

  const bool rm2xx2iszero      = rm2xx2 == ScalarT(0.0);
  const bool rm2xy_rm2xyiszero = rm2xy_rm2xy == ScalarT(0.0);
  const bool both_zero         = rm2xx2iszero && rm2xy_rm2xyiszero;

  // check degeneracy
  evec0[0] = if_then_else(both_zero, a_row2[0], evec0[0]);
  evec0[1] = if_then_else(both_zero, a_row2[1], evec0[1]);
  evec0[2] = if_then_else(both_zero, a_row2[2], evec0[2]);

  evec1[K_X] = evec2[K_Y] * evec0[K_Z] - evec2[K_Z] * evec0[K_Y];
  evec1[K_Y] = evec2[K_Z] * evec0[K_X] - evec2[K_X] * evec0[K_Z];
  evec1[K_Z] = evec2[K_X] * evec0[K_Y] - evec2[K_Y] * evec0[K_X];

  eval0 += c1;
  eval1 += c1;
  eval2 += c1;

  return;
}

template <typename ScalarT>
NIMBLE_INLINE_FUNCTION void
Polar_Decomp(
    const ScalarT* const mat,          /* full tensor */
    ScalarT* const       left_stretch, /* symmetric tensor */
    ScalarT* const       rotation)
{ /* full tensor */

  ScalarT mat_inv[9];
  Invert_Full33(mat, mat_inv);

  ScalarT mat_inv_squared[6];
  Square_Full33T_Full33(mat_inv, mat_inv_squared);

  ScalarT vec1[3] = {0.0, 0.0, 0.0};  // Eigen vectors of mat_inv_squared
  ScalarT vec2[3] = {0.0, 0.0, 0.0};
  ScalarT vec3[3] = {0.0, 0.0, 0.0};
  ScalarT eval[3] = {0.0, 0.0, 0.0};  // Eigen values of mat_inv_squared

  // Eigen_Sym33_NonUnit(mat_inv_squared, eval[0], eval[1], eval[2], vec1, vec2, vec3);
  Eigen_Sym33(mat_inv_squared, eval[0], eval[1], eval[2], vec1, vec2, vec3);

  const ScalarT zero = ScalarT(0.0);

  eval[0] = if_then(eval[0] < zero, zero, eval[0]);
  eval[1] = if_then(eval[1] < zero, zero, eval[1]);
  eval[2] = if_then(eval[2] < zero, zero, eval[2]);

  const ScalarT len1_sq = vec1[K_X] * vec1[K_X] + vec1[K_Y] * vec1[K_Y] + vec1[K_Z] * vec1[K_Z];
  const ScalarT len2_sq = vec2[K_X] * vec2[K_X] + vec2[K_Y] * vec2[K_Y] + vec2[K_Z] * vec2[K_Z];
  const ScalarT len3_sq = vec3[K_X] * vec3[K_X] + vec3[K_Y] * vec3[K_Y] + vec3[K_Z] * vec3[K_Z];

  const ScalarT xlx = std::sqrt(eval[0]);
  const ScalarT xly = std::sqrt(eval[1]);
  const ScalarT xlz = std::sqrt(eval[2]);

  const ScalarT one = ScalarT(1.0);

  const ScalarT xlxi = one / (xlx * len1_sq);
  const ScalarT xlyi = one / (xly * len2_sq);
  const ScalarT xlzi = one / (xlz * len3_sq);

  left_stretch[K_S_XX] = xlxi * vec1[K_X] * vec1[K_X] + xlyi * vec2[K_X] * vec2[K_X] + xlzi * vec3[K_X] * vec3[K_X];
  left_stretch[K_S_YY] = xlxi * vec1[K_Y] * vec1[K_Y] + xlyi * vec2[K_Y] * vec2[K_Y] + xlzi * vec3[K_Y] * vec3[K_Y];
  left_stretch[K_S_ZZ] = xlxi * vec1[K_Z] * vec1[K_Z] + xlyi * vec2[K_Z] * vec2[K_Z] + xlzi * vec3[K_Z] * vec3[K_Z];
  left_stretch[K_S_XY] = xlxi * vec1[K_X] * vec1[K_Y] + xlyi * vec2[K_X] * vec2[K_Y] + xlzi * vec3[K_X] * vec3[K_Y];
  left_stretch[K_S_YZ] = xlxi * vec1[K_Y] * vec1[K_Z] + xlyi * vec2[K_Y] * vec2[K_Z] + xlzi * vec3[K_Y] * vec3[K_Z];
  left_stretch[K_S_ZX] = xlxi * vec1[K_Z] * vec1[K_X] + xlyi * vec2[K_Z] * vec2[K_X] + xlzi * vec3[K_Z] * vec3[K_X];

  Mult_Full33_Sym33_ReturnT(mat_inv, left_stretch, rotation);
}

template <typename ScalarT>
NIMBLE_INLINE_FUNCTION void
Polar_Left_LogV_Lame(
    const ScalarT* const def_grad_inc,     /* full tensor */
    ScalarT* const       left_stretch_inc, /* symmetric tensor */
    ScalarT* const       rotation_inc,     /* full tensor */
    ScalarT* const       log_left_stretch_inc)
{ /* symmetric tensor */

  // Adapted from MiniTensor_LinearAlgebra.t.h, polar_left_logV_lame()

  // compute spd tensor
  ScalarT b[6];
  Square_Full33_Full33T(def_grad_inc, b);

  // get eigenvalues/eigenvectors
  ScalarT vec1[3] = {0.0, 0.0, 0.0};  // Eigen vectors of b
  ScalarT vec2[3] = {0.0, 0.0, 0.0};
  ScalarT vec3[3] = {0.0, 0.0, 0.0};
  ScalarT eval[3] = {0.0, 0.0, 0.0};  // Eigen values of b
  // Eigen_Sym33_NonUnit(b, eval[0], eval[1], eval[2], vec1, vec2, vec3);
  Eigen_Sym33(b, eval[0], eval[1], eval[2], vec1, vec2, vec3);

  CheckVectorSanity(3, eval, "Eigenvalues in Polar_Left_LogV_Lame()");

  // normalize the eigenvectors
  ScalarT mag;
  mag = std::sqrt(vec1[0] * vec1[0] + vec1[1] * vec1[1] + vec1[2] * vec1[2]);
  vec1[0] /= mag;
  vec1[1] /= mag;
  vec1[2] /= mag;
  mag = std::sqrt(vec2[0] * vec2[0] + vec2[1] * vec2[1] + vec2[2] * vec2[2]);
  vec2[0] /= mag;
  vec2[1] /= mag;
  vec2[2] /= mag;
  mag = std::sqrt(vec3[0] * vec3[0] + vec3[1] * vec3[1] + vec3[2] * vec3[2]);
  vec3[0] /= mag;
  vec3[1] /= mag;
  vec3[2] /= mag;

  // std::cout << "DJL DEBUGGING eigen values in Polar_Left_LogV_Lame " <<
  // eval[0] << ", " << eval[1] << ", " << eval[2] << std::endl; std::cout <<
  // "DJL DEBUGGING eigen vectors in Polor_Left_LogV_Lame (" << vec1[0] << ", "
  // << vec1[1] << ", " << vec1[2] << ") (" << vec2[0] << ", " << vec2[1] << ",
  // "
  // << vec2[2] << ") (" << vec3[0] << ", " << vec3[1] << ", " << vec3[2] << ")"
  // << std::endl;

  ScalarT x[3];
  x[0] = std::sqrt(eval[0]);
  x[1] = std::sqrt(eval[1]);
  x[2] = std::sqrt(eval[2]);

  ScalarT xi[3];
  xi[0] = 1.0 / x[0];
  xi[1] = 1.0 / x[1];
  xi[2] = 1.0 / x[2];

  ScalarT lnx[3];
  lnx[0] = std::log(x[0]);
  lnx[1] = std::log(x[1]);
  lnx[2] = std::log(x[2]);

  ScalarT vec_transpose[9];
  vec_transpose[K_F_XX] = vec1[0];
  vec_transpose[K_F_XY] = vec1[1];
  vec_transpose[K_F_XZ] = vec1[2];
  vec_transpose[K_F_YX] = vec2[0];
  vec_transpose[K_F_YY] = vec2[1];
  vec_transpose[K_F_YZ] = vec2[2];
  vec_transpose[K_F_ZX] = vec3[0];
  vec_transpose[K_F_ZY] = vec3[1];
  vec_transpose[K_F_ZZ] = vec3[2];

  ScalarT temp[9];

  // V = eVec*x*transpose(eVec)
  temp[K_F_XX]             = vec1[0] * x[0];
  temp[K_F_XY]             = vec2[0] * x[1];
  temp[K_F_XZ]             = vec3[0] * x[2];
  temp[K_F_YX]             = vec1[1] * x[0];
  temp[K_F_YY]             = vec2[1] * x[1];
  temp[K_F_YZ]             = vec3[1] * x[2];
  temp[K_F_ZX]             = vec1[2] * x[0];
  temp[K_F_ZY]             = vec2[2] * x[1];
  temp[K_F_ZZ]             = vec3[2] * x[2];
  left_stretch_inc[K_S_XX] = temp[K_F_XX] * vec_transpose[K_F_XX] + temp[K_F_XY] * vec_transpose[K_F_YX] +
                             temp[K_F_XZ] * vec_transpose[K_F_ZX];
  left_stretch_inc[K_S_YY] = temp[K_F_YX] * vec_transpose[K_F_XY] + temp[K_F_YY] * vec_transpose[K_F_YY] +
                             temp[K_F_YZ] * vec_transpose[K_F_ZY];
  left_stretch_inc[K_S_ZZ] = temp[K_F_ZX] * vec_transpose[K_F_XZ] + temp[K_F_ZY] * vec_transpose[K_F_YZ] +
                             temp[K_F_ZZ] * vec_transpose[K_F_ZZ];
  left_stretch_inc[K_S_XY] = temp[K_F_XX] * vec_transpose[K_F_XY] + temp[K_F_XY] * vec_transpose[K_F_YY] +
                             temp[K_F_XZ] * vec_transpose[K_F_ZY];
  left_stretch_inc[K_S_YZ] = temp[K_F_YX] * vec_transpose[K_F_XZ] + temp[K_F_YY] * vec_transpose[K_F_YZ] +
                             temp[K_F_YZ] * vec_transpose[K_F_ZZ];
  left_stretch_inc[K_S_ZX] = temp[K_F_ZX] * vec_transpose[K_F_XX] + temp[K_F_ZY] * vec_transpose[K_F_YX] +
                             temp[K_F_ZZ] * vec_transpose[K_F_ZX];

  // Vinv = eVec*xi*transpose(eVec)
  temp[K_F_XX] = vec1[0] * xi[0];
  temp[K_F_XY] = vec2[0] * xi[1];
  temp[K_F_XZ] = vec3[0] * xi[2];
  temp[K_F_YX] = vec1[1] * xi[0];
  temp[K_F_YY] = vec2[1] * xi[1];
  temp[K_F_YZ] = vec3[1] * xi[2];
  temp[K_F_ZX] = vec1[2] * xi[0];
  temp[K_F_ZY] = vec2[2] * xi[1];
  temp[K_F_ZZ] = vec3[2] * xi[2];
  ScalarT left_stretch_inc_inverse[9];
  Mult_Full33_Full33(temp, vec_transpose, left_stretch_inc_inverse);

  // Vinv = eVec*lnx*transpose(eVec)
  temp[K_F_XX]                 = vec1[0] * lnx[0];
  temp[K_F_XY]                 = vec2[0] * lnx[1];
  temp[K_F_XZ]                 = vec3[0] * lnx[2];
  temp[K_F_YX]                 = vec1[1] * lnx[0];
  temp[K_F_YY]                 = vec2[1] * lnx[1];
  temp[K_F_YZ]                 = vec3[1] * lnx[2];
  temp[K_F_ZX]                 = vec1[2] * lnx[0];
  temp[K_F_ZY]                 = vec2[2] * lnx[1];
  temp[K_F_ZZ]                 = vec3[2] * lnx[2];
  log_left_stretch_inc[K_S_XX] = temp[K_F_XX] * vec_transpose[K_F_XX] + temp[K_F_XY] * vec_transpose[K_F_YX] +
                                 temp[K_F_XZ] * vec_transpose[K_F_ZX];
  log_left_stretch_inc[K_S_YY] = temp[K_F_YX] * vec_transpose[K_F_XY] + temp[K_F_YY] * vec_transpose[K_F_YY] +
                                 temp[K_F_YZ] * vec_transpose[K_F_ZY];
  log_left_stretch_inc[K_S_ZZ] = temp[K_F_ZX] * vec_transpose[K_F_XZ] + temp[K_F_ZY] * vec_transpose[K_F_YZ] +
                                 temp[K_F_ZZ] * vec_transpose[K_F_ZZ];
  log_left_stretch_inc[K_S_XY] = temp[K_F_XX] * vec_transpose[K_F_XY] + temp[K_F_XY] * vec_transpose[K_F_YY] +
                                 temp[K_F_XZ] * vec_transpose[K_F_ZY];
  log_left_stretch_inc[K_S_YZ] = temp[K_F_YX] * vec_transpose[K_F_XZ] + temp[K_F_YY] * vec_transpose[K_F_YZ] +
                                 temp[K_F_YZ] * vec_transpose[K_F_ZZ];
  log_left_stretch_inc[K_S_ZX] = temp[K_F_ZX] * vec_transpose[K_F_XX] + temp[K_F_ZY] * vec_transpose[K_F_YX] +
                                 temp[K_F_ZZ] * vec_transpose[K_F_ZX];

  // R = Vinv*F
  Mult_Full33_Full33(left_stretch_inc_inverse, def_grad_inc, rotation_inc);
}

template <typename ScalarT>
NIMBLE_INLINE_FUNCTION void
CheckCorrectnessOfPolarDecomp(
    const ScalarT* const mat,
    const ScalarT* const V,
    const ScalarT* const R,
    std::string const&   label)
{
  ScalarT norm_tol = 1.0e-12;

  ScalarT check[9], should_be_zero[9];
  Mult_Sym33_Full33(V, R, check);
  for (int i = 0; i < 9; i++) { should_be_zero[i] = mat[i] - check[i]; }
  ScalarT should_be_zero_norm = Norm_Full33(should_be_zero);
  if (should_be_zero_norm > norm_tol) {
    std::cout << "\n\nFailure in CheckCorrectnessOfPolarDecomp()!" << std::endl;
    Print_Full33("mat", mat);
    Print_Full33("V", V);
    Print_Full33("R", R);
    std::cout << "norm of error = " << should_be_zero_norm << "\n" << std::endl;
    // NIMBLE_ABORT("\n**** Polar decomposition correctness check
    // failed for " + label + "!\n");
  }
}

template <typename ScalarT>
NIMBLE_INLINE_FUNCTION void
CheckCorrectnessOfSymSkew(
    const ScalarT* const mat,
    const ScalarT* const sym,
    const ScalarT* const skew,
    std::string const&   label)
{
  ScalarT norm_tol = 1.0e-12 * Norm_Full33(mat);

  ScalarT should_be_zero[9];
  should_be_zero[K_F_XX]      = sym[K_S_XX] + skew[K_F_XX] - mat[K_F_XX];
  should_be_zero[K_F_XY]      = sym[K_S_XY] + skew[K_F_XY] - mat[K_F_XY];
  should_be_zero[K_F_XZ]      = sym[K_S_XZ] + skew[K_F_XZ] - mat[K_F_XZ];
  should_be_zero[K_F_YX]      = sym[K_S_YX] + skew[K_F_YX] - mat[K_F_YX];
  should_be_zero[K_F_YY]      = sym[K_S_YY] + skew[K_F_YY] - mat[K_F_YY];
  should_be_zero[K_F_YZ]      = sym[K_S_YZ] + skew[K_F_YZ] - mat[K_F_YZ];
  should_be_zero[K_F_ZX]      = sym[K_S_ZX] + skew[K_F_ZX] - mat[K_F_ZX];
  should_be_zero[K_F_ZY]      = sym[K_S_ZY] + skew[K_F_ZY] - mat[K_F_ZY];
  should_be_zero[K_F_ZZ]      = sym[K_S_ZZ] + skew[K_F_ZZ] - mat[K_F_ZZ];
  ScalarT should_be_zero_norm = Norm_Full33(should_be_zero);
  if (should_be_zero_norm > norm_tol) {
    std::cout << "\n\nFailure in CheckCorrectnessOfSymSkew()!" << std::endl;
    Print_Full33("mat", mat);
    Print_Sym33("sym", sym);
    Print_Full33("skew", skew);
    std::cout << "norm of error = " << should_be_zero_norm << "\n" << std::endl;
    // NIMBLE_ABORT("\n**** Sym-skew correctness check failed for " +
    // label + "!\n");
  }
}

template <typename ScalarT>
NIMBLE_INLINE_FUNCTION void
Log_Rotation_Pi(const ScalarT* const rotation, ScalarT* const log_rotation)
{
  // adapted from MiniTensor_LinearAlgebra.t.h log_rotation_pi()

  ScalarT machine_epsilon = std::numeric_limits<ScalarT>::epsilon();

  // set firewall to make sure the rotation is indeed 180 degrees
  if (std::abs(rotation[K_F_XX] + rotation[K_F_YY] + rotation[K_F_ZZ] + 1.0) < 10.0 * machine_epsilon) {
    // throw std::invalid_argument("\n**** Input check failed for Log_Rotation_Pi()!\n");
  }

  ScalarT B[9];
  B[K_F_XX] = rotation[K_F_XX] - 1.0;
  B[K_F_XY] = rotation[K_F_XY];
  B[K_F_XZ] = rotation[K_F_XZ];
  B[K_F_YX] = rotation[K_F_YX];
  B[K_F_YY] = rotation[K_F_YY] - 1.0;
  B[K_F_YZ] = rotation[K_F_YZ];
  B[K_F_ZX] = rotation[K_F_ZX];
  B[K_F_ZY] = rotation[K_F_ZY];
  B[K_F_ZZ] = rotation[K_F_ZZ] - 1.0;

  ScalarT u[3];
  u[0] = B[K_F_XX];
  u[1] = B[K_F_XY];
  u[2] = B[K_F_XZ];

  ScalarT v[3];
  v[0] = B[K_F_YX];
  v[1] = B[K_F_YY];
  v[2] = B[K_F_YZ];

  ScalarT normal[3];
  CrossProduct(u, v, normal);

  ScalarT norm = normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2];
  if (norm > 0.0) { norm = std::sqrt(norm); }

  if (norm < machine_epsilon) {
    ScalarT w[3];
    w[0] = B[K_F_ZX];
    w[1] = B[K_F_ZY];
    w[2] = B[K_F_ZZ];

    CrossProduct(u, w, normal);

    norm = normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2];
    if (norm > 0.0) { norm = std::sqrt(norm); }

    NIMBLE_ASSERT(norm >= machine_epsilon, "\n**** Error in Log_Rotation_Pi(), cannot"
                  " determine rotation vector.");
  }

  normal[0] /= norm;
  normal[1] /= norm;
  normal[2] /= norm;

  ScalarT pi = std::acos(-1.0);

  log_rotation[K_F_XX] = 0.0;
  log_rotation[K_F_XY] = -normal[2] * pi;
  log_rotation[K_F_XZ] = normal[1] * pi;
  log_rotation[K_F_YX] = normal[2] * pi;
  log_rotation[K_F_YY] = 0.0;
  log_rotation[K_F_YZ] = -normal[0] * pi;
  log_rotation[K_F_ZX] = -normal[1] * pi;
  log_rotation[K_F_ZY] = normal[0] * pi;
  log_rotation[K_F_ZZ] = 0.0;
}

template <typename ScalarT>
NIMBLE_INLINE_FUNCTION void
Log_Rotation(const ScalarT* const rotation, ScalarT* const log_rotation)
{
  // adapted from MiniTensor_LinearAlgebra log_rotation()

  // firewalls, make sure R \in SO(N)
  bool    valid_input     = true;
  ScalarT machine_epsilon = std::numeric_limits<ScalarT>::epsilon();
  ScalarT temp[6];
  Square_Full33_Full33T(rotation, temp);
  temp[K_S_XX] -= 1.0;
  temp[K_S_YY] -= 1.0;
  temp[K_S_ZZ] -= 1.0;
  ScalarT norm = Norm_Sym33(temp);
  if (norm > 100.0 * machine_epsilon) { valid_input = false; }
  ScalarT det_rotation = Determinant_Full33(rotation);
  if (std::abs(det_rotation - 1.0) > 100.0 * machine_epsilon) { valid_input = false; }
  if (!valid_input) {
    std::stringstream ss;
    ss << "\n**** Input check failed for Log_Rotation()!";
    ss << "\n****   norm(R - I)  must be less than (100 * machine_epsilon):  "
          "norm = "
       << norm << ", machine_epsilon = " << machine_epsilon;
    ss << "\n****   (det(R) - 1) must be less than (100 * machine_epsilon):  "
          "det(R) = "
       << det_rotation << ", machine_epsilon = " << machine_epsilon << "\n";
    // throw std::invalid_argument(ss.str());
  }

  // acos requires input between -1 and +1
  ScalarT cosine = 0.5 * (rotation[K_F_XX] + rotation[K_F_YY] + rotation[K_F_ZZ] - 1.0);
  if (cosine < -1.0) {
    cosine = -1.0;
  } else if (cosine > 1.0) {
    cosine = 1.0;
  }
  ScalarT theta = std::acos(cosine);

  if (theta == 0.0) {
    Zero_Full33(log_rotation);
  } else if (std::abs(cosine + 1.0) < 10.0 * machine_epsilon) {
    //        r = log_rotation_pi(R);
  } else {
    // r = theta / std::sin(theta) * skew(R);
    ScalarT a = theta / std::sin(theta);
    ScalarT rotation_skew[9];
    SkewPart_Full33(rotation, rotation_skew);
    Mult_Scalar_Full33(a, rotation_skew, log_rotation);
  }
}

NIMBLE_INLINE_FUNCTION double
Invert3x3(const double mat[][3], double inv[][3])
{
  double minor0 = mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1];
  double minor1 = mat[1][0] * mat[2][2] - mat[1][2] * mat[2][0];
  double minor2 = mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0];
  double minor3 = mat[0][1] * mat[2][2] - mat[0][2] * mat[2][1];
  double minor4 = mat[0][0] * mat[2][2] - mat[2][0] * mat[0][2];
  double minor5 = mat[0][0] * mat[2][1] - mat[0][1] * mat[2][0];
  double minor6 = mat[0][1] * mat[1][2] - mat[0][2] * mat[1][1];
  double minor7 = mat[0][0] * mat[1][2] - mat[0][2] * mat[1][0];
  double minor8 = mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0];
  double det    = mat[0][0] * minor0 - mat[0][1] * minor1 + mat[0][2] * minor2;

  if (det <= 0.0) {
#ifdef NIMBLE_HAVE_KOKKOS
    if (det == 0.0)
      printf("\n**** Error in Invert3x3(), singular matrix.\n");
    else
      printf(
          "\n**** Error in HexInvert3x3(), negative determinant "
          "(%e)\n",
          det);
#else
    NIMBLE_ASSERT(det > 0.0, "\n**** Error in Invert3x3(), singular matrix.\n");
#endif
  }

  inv[0][0] = minor0 / det;
  inv[0][1] = -1.0 * minor3 / det;
  inv[0][2] = minor6 / det;
  inv[1][0] = -1.0 * minor1 / det;
  inv[1][1] = minor4 / det;
  inv[1][2] = -1.0 * minor7 / det;
  inv[2][0] = minor2 / det;
  inv[2][1] = -1.0 * minor5 / det;
  inv[2][2] = minor8 / det;

  return det;
}


template< class Matrix >
void
LU_Decompose(int num_entries, Matrix& mat, int* index)
{
  int                 imax = 0;
  double              big, temp, sum;
  int                 n = num_entries;
  std::vector<double> vv(n);  // vv stores the implicit scaling of each row.
  double              tiny = 1e-40;

  // loop over rows to get the implicit scaling information
  // scale each row by 1/(largest element)

  for (int i = 0; i < n; ++i) {
    big = 0.0;
    for (int j = 0; j < n; ++j) {
      if (std::fabs(mat(i, j)) > big) big = std::fabs(mat(i, j));
    }
    if (big < tiny) {
      // No nonzero largest element.
      NIMBLE_ABORT("Singular matrix in routine LU_Decompose()");
    }
    vv[i] = 1.0 / big;
  }

  // loop over columns
  for (int j = 0; j < n; ++j) {
    for (int i = 0; i < j; ++i) {
      sum = mat(i, j);
      for (int k = 0; k < i; ++k) { sum -= mat(i, k) * mat(k, j); }
      mat(i, j) = sum;
    }

    big = 0.0;
    for (int i = j; i < n; ++i) {
      sum = mat(i, j);
      for (int k = 0; k < j; ++k) { sum -= mat(i, k) * mat(k, j); }
      mat(i, j) = sum;

      if (vv[i] * std::fabs(sum) > big) {
        // this is the biggest element in the column
        // save it as the pivot element
        big  = vv[i] * std::fabs(sum);
        imax = i;
      }
    }

    // interchange rows if required
    if (j != imax) {
      for (int k = 0; k < n; ++k) {
        temp         = mat(imax, k);
        mat(imax, k) = mat(j, k);
        mat(j, k)    = temp;
      }
      vv[imax] = vv[j];
    }
    index[j] = imax;

    if (std::fabs(mat(j, j)) < tiny) {
      // matrix is singular and we're about to divide by zero
      NIMBLE_ABORT("Singular matrix in routine LU_Decompose()");
    }

    // divide by the pivot element
    if (j != n) {
      temp = 1.0 / (mat(j, j));
      for (int i = j + 1; i < n; ++i) { mat(i, j) *= temp; }
    }
  }

  // decomposition is returned in mat
}


template< class Matrix >
void
LU_Solve(int num_entries, const Matrix& mat, double* vec, const int* index)
{
  // mat holds the upper and lower triangular decomposition (obtained by
  // LU_Deompose) vec hold the right-hand-side vector index holds the pivoting
  // information (row interchanges)

  double sum;
  int    n = num_entries;

  // forward substitution
  // resolve pivoting as we go (stored in index array)
  for (int i = 0; i < n; ++i) {
    sum           = vec[index[i]];
    vec[index[i]] = vec[i];
    for (int j = 0; j < i; ++j) { sum -= mat(i, j) * vec[j]; }
    vec[i] = sum;
  }

  // back substitution
  for (int i = n - 1; i >= 0; --i) {
    sum = vec[i];
    for (int j = i + 1; j < n; ++j) { sum -= mat(i, j) * vec[j]; }
    vec[i] = sum / mat(i, i);
  }

  // solution is returned in vec
}


template< class Matrix >
void
LU_SolveSystem(int num_entries, Matrix& mat, double* vec, int* scratch)
{
  int* index = scratch;

  // decompose the matrix into upper and lower triangular parts
  LU_Decompose(num_entries, mat, index);

  // solve the system (forward and backward substitution)
  LU_Solve(num_entries, mat, vec, index);
}

#endif  // NIMBLE_UTILS_H
