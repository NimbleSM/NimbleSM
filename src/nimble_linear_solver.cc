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

#include "nimble_linear_solver.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "nimble_macros.h"

namespace nimble {

void
CRSMatrixContainer::AllocateNonzeros(std::vector<int> const& i_index, std::vector<int> const& j_index)
{
  i_index_ = i_index;
  j_index_ = j_index;
  data_.resize(i_index.size(), 0.0);
  int max_row = 0;
  for (auto const& i : i_index) {
    if (i > max_row) { max_row = i; }
  }
  num_rows_ = max_row + 1;
  row_first_index_.resize(num_rows_ + 1);
  int current_row = -1;
  for (int i = 0; i < i_index_.size(); i++) {
    int row = i_index_[i];
    if (row > current_row) {
      row_first_index_[row] = i;
      current_row           = row;
    }
  }
  row_first_index_[num_rows_] = i_index_.size();
}

void
CRSMatrixContainer::SetRowValues(int row, double value)
{
  for (int i = row_first_index_[row]; i < row_first_index_[row + 1]; i++) { data_[i] = value; }
}

void
CRSMatrixContainer::SetColumnValues(int col, double value)
{
  for (int i = 0; i < j_index_.size(); i++) {
    if (j_index_[i] == col) { data_[i] = value; }
  }
}

int
CRSMatrixContainer::FindIndex(int i_index, int j_index) const
{
  int min        = row_first_index_[i_index];
  int max        = row_first_index_[i_index + 1];
  int trial      = (min + max) / 2;
  int j_value    = j_index_[trial];
  int max_trials = max - min;

  for (int i = 0; i < max_trials; i++) {
    if (j_value == j_index) {
      return trial;
    } else if (j_value < j_index) {
      min = trial;
    } else {
      max = trial;
    }
    trial   = (min + max) / 2;
    j_value = j_index_[trial];
  }

  throw std::invalid_argument("Error, CRSMatrixContainer::FindIndex() failed to find index.");
}

void
CRSMatrixContainer::MatVec(const double* vec, double* result) const
{
  for (unsigned int i_row = 0; i_row < num_rows_; i_row++) { result[i_row] = 0.0; }
  for (unsigned int i = 0; i < data_.size(); i++) { result[i_index_[i]] += data_[i] * vec[j_index_[i]]; }
}

void
CRSMatrixContainer::DiagonalMatrixMatVec(const double* vec, double* result) const
{
  // the matrix must be diagonal
  if (num_rows_ != data_.size()) {
    NIMBLE_ABORT(
        "**** Error in CRSMatrixContainer::DiagonalMatrixMatVec(), matrix is "
        "not diagonal.\n");
  }

  for (unsigned int i = 0; i < num_rows_; i++) { result[i] = data_[i] * vec[i]; }
}

void
PopulateDiagonalPreconditioner(const CRSMatrixContainer& A, CRSMatrixContainer& M)
{
  // the preconditioner must be a diagonal matrix
  if (M.num_rows_ != M.data_.size()) {
    throw std::invalid_argument(
        "**** Error in PopulateDiagonalPreconditioner(), preconditioner must "
        "be a diagonal matrix.\n");
  }

  int m_index(0);
  for (int i = 0; i < A.data_.size(); i++) {
    if (A.i_index_[i] == A.j_index_[i]) { M.data_[m_index++] = 1.0 / A.data_[i]; }
  }
}


bool
CG_SolveSystem(
    nimble::CRSMatrixContainer& A,
    const double*               b,
    CGScratchSpace&             cg_scratch,
    double*                     x,
    int&                        num_iterations,
    double                      cg_tol,
    int                         max_iterations)
{
  // see "An Introduction to the Conjugate Gradient Method Without the Agonizing
  // Pain", J.R. Shewchuk, 1994.

  double       alpha, beta, delta_old, delta_new;
  unsigned int num_entries = A.NumRows();
  cg_scratch.Resize(num_entries);
  double* d = cg_scratch.d.data();
  double* r = cg_scratch.r.data();
  double* s = cg_scratch.s.data();
  double* q = cg_scratch.q.data();

  // diagonal preconditioner
  CRSMatrixContainer& M = cg_scratch.M;
  PopulateDiagonalPreconditioner(A, M);

  // r = b - Ax
  A.MatVec(x, q);
  for (unsigned int i = 0; i < num_entries; i++) { r[i] = b[i] - q[i]; }

  // d = M^-1 r
  M.DiagonalMatrixMatVec(r, d);

  // delta_new = r^T d
  // delta_old = delta_new
  delta_new = delta_old = InnerProduct(num_entries, r, d);
  double tolerance      = cg_tol * delta_old;

  int iteration      = 0;

  while (delta_new > tolerance && iteration < max_iterations) {
    // q = Ad
    A.MatVec(d, q);

    // alpha = delta_new / (d^T q)
    alpha = delta_new / InnerProduct(num_entries, d, q);

    // x = x + alpha * d
    for (unsigned int i = 0; i < num_entries; i++) { x[i] += alpha * d[i]; }

    if (iteration % 50 == 0) {
      // r = b - Ax
      A.MatVec(x, q);  // here, q is just a place to store Ax
      for (unsigned int i = 0; i < num_entries; i++) { r[i] = b[i] - q[i]; }
    } else {
      // r = r - alpha * q
      for (unsigned int i = 0; i < num_entries; i++) { r[i] -= alpha * q[i]; }
    }

    // s = M^-1 r
    M.DiagonalMatrixMatVec(r, s);

    // delta_old = delta_new
    delta_old = delta_new;

    // delta_new = r^T s
    delta_new = InnerProduct(num_entries, r, s);

    // beta = delta_new / delta_old
    beta = delta_new / delta_old;

    // d = s + beta * d
    for (unsigned int i = 0; i < num_entries; i++) { d[i] = s[i] + beta * d[i]; }

    iteration += 1;
  }

  num_iterations = iteration;
  bool success   = true;
  if (iteration == max_iterations) { success = false; }
  return success;
}

}  // namespace nimble
