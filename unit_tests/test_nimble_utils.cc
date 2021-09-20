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

#include <vector>

#include <gtest/gtest.h>
#include "nimble_utils.h"

TEST(nimble_utils, invert_3x3_matrix)
{

  std::vector<double> mat(9);
  std::vector<double> inv(9), ref_inv(9);

  mat[K_F_XX] = 1.0; mat[K_F_XY] = 0.0; mat[K_F_XZ] = 0.0;
  mat[K_F_YX] = 1.0; mat[K_F_YY] = 1.0; mat[K_F_YZ] = 0.0;
  mat[K_F_ZX] = 1.0; mat[K_F_ZY] = 1.0; mat[K_F_ZZ] = 1.0;

  ref_inv[K_F_XX] = 1.0; ref_inv[K_F_XY] = 0.0; ref_inv[K_F_XZ] = 0.0;
  ref_inv[K_F_YX] = -1.0; ref_inv[K_F_YY] = 1.0; ref_inv[K_F_YZ] = 0.0;
  ref_inv[K_F_ZX] = 0.0; ref_inv[K_F_ZY] = -1.0; ref_inv[K_F_ZZ] = 1.0;

  double det = Invert_Full33(&mat[0], &inv[0]);

  bool correctInverse = true;
  for (int ii = 0; ii < 9; ++ii) {
    if (ref_inv[ii] != inv[ii]) {
      correctInverse = false;
      break;
    }
  }

  ASSERT_TRUE(correctInverse);
  ASSERT_DOUBLE_EQ(det, 1.0);

  mat[K_F_XX] = 3.0; mat[K_F_XY] = 0.0; mat[K_F_XZ] = 2.0;
  mat[K_F_YX] = 2.0; mat[K_F_YY] = 0.0; mat[K_F_YZ] = -2.0;
  mat[K_F_ZX] = 0.0; mat[K_F_ZY] = 1.0; mat[K_F_ZZ] = 1.0;

  ref_inv[K_F_XX] = 0.2; ref_inv[K_F_XY] = 0.2; ref_inv[K_F_XZ] = 0.0;
  ref_inv[K_F_YX] = -0.2; ref_inv[K_F_YY] = 0.3; ref_inv[K_F_YZ] = 1.0;
  ref_inv[K_F_ZX] = 0.2; ref_inv[K_F_ZY] = -0.3; ref_inv[K_F_ZZ] = 0.0;

  det = Invert_Full33(&mat[0], &inv[0]);

  correctInverse = true;
  for (int ii = 0; ii < 9; ++ii) {
    if (ref_inv[ii] != inv[ii]) {
      correctInverse = false;
      break;
    }
  }

  ASSERT_TRUE(correctInverse);
  ASSERT_DOUBLE_EQ(det, 10.0);

}
