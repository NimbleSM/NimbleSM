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

#ifndef NIMBLE_LINEAR_SOLVER_H
#define NIMBLE_LINEAR_SOLVER_H

#ifdef NIMBLE_HAVE_DARMA
#include "darma.h"
#else
#include <vector>
#endif

#ifdef NIMBLE_HAVE_MPI
#include <mpi.h>
#endif

namespace nimble {

class MatrixContainer
{
 public:
  MatrixContainer() : num_rows_(0) {}

  ~MatrixContainer() {}

  void
  AllocateNonzeros(std::vector<int> const& i_index, std::vector<int> const& j_index)
  {
    num_rows_ = 0;
    for (auto const& entry : i_index) {
      if (entry > num_rows_) { num_rows_ = entry; }
    }
    num_rows_ += 1;
    data_.resize(num_rows_ * num_rows_, 0.0);
  }

  void
  SetAllValues(double value)
  {
    for (unsigned int i = 0; i < data_.size(); i++) { data_[i] = value; }
  }

  void
  SetRowValues(int row, double value)
  {
    for (int i = 0; i < num_rows_; ++i) { data_[num_rows_ * row + i] = value; }
  }

  void
  SetColumnValues(int col, double value)
  {
    for (int i = 0; i < num_rows_; ++i) { data_[num_rows_ * i + col] = value; }
  }

  unsigned int
  NumNonzeros() const
  {
    return data_.size();
  }

  inline double
  operator()(int i, int j) const
  {
    return data_[num_rows_ * i + j];
  }

  inline double&
  operator()(int i, int j)
  {
    return data_[num_rows_ * i + j];
  }

 private:
  int                 num_rows_;
  std::vector<double> data_;
};

class CRSMatrixContainer
{
 public:
  CRSMatrixContainer() : num_rows_(0) {}

  ~CRSMatrixContainer() = default;

  void
  AllocateNonzeros(std::vector<int> const& i_index, std::vector<int> const& j_index);

  void
  AllocateDiagonalMatrix(int num_rows)
  {
    if (num_rows_ == num_rows) { return; }
    num_rows_ = num_rows;
    data_.resize(num_rows_);
    i_index_.resize(num_rows_);
    j_index_.resize(num_rows_);
    row_first_index_.resize(num_rows_);
    for (int i = 0; i < num_rows_; i++) {
      data_[i]            = 0.0;
      i_index_[i]         = i;
      j_index_[i]         = i;
      row_first_index_[i] = i;
    }
  }

  void
  SetAllValues(double value)
  {
    for (double& d_data : data_) d_data = value;
  }

  void
  SetRowValues(int row, double value);

  void
  SetColumnValues(int col, double value);

  int
  NumRows() const
  {
    return num_rows_;
  }

  unsigned int
  NumNonzeros() const
  {
    return data_.size();
  }

  void
  MatVec(const double* vec, double* result) const;

  void
  DiagonalMatrixMatVec(const double* vec, double* result) const;

  inline double
  operator()(int i_index, int j_index) const
  {
    return data_[FindIndex(i_index, j_index)];
  }

  inline double&
  operator()(int i_index, int j_index)
  {
    return data_[FindIndex(i_index, j_index)];
  }

#ifdef NIMBLE_HAVE_DARMA
  template <typename ArchiveType>
  void
  serialize(ArchiveType& ar)
  {
    ar | num_rows_ | data_ | i_index_ | j_index_ | row_first_index_;
  }
#endif

 private:
  int
  FindIndex(int i_index, int j_index) const;

  int                 num_rows_;
  std::vector<double> data_;
  std::vector<int>    i_index_;
  std::vector<int>    j_index_;
  std::vector<int>    row_first_index_;

  /// \brief Generate a diagonal preconditioner from a sparse matrix
  ///
  /// \param A Input sparse matrix
  /// \param M Diagonal preconditioner (mathematically M = (diag(A))^{-1})
  friend void
  PopulateDiagonalPreconditioner(const CRSMatrixContainer& A, CRSMatrixContainer& M);
};

inline double
InnerProduct(unsigned int num_entries, const double* vec_1, const double* vec_2)
{
  double result(0.0);
  for (unsigned int i = 0; i < num_entries; i++) { result += vec_1[i] * vec_2[i]; }
#ifdef NIMBLE_HAVE_MPI
  double restmp = result;
  MPI_Allreduce(&restmp, &result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
  return result;
}

inline double
InnerProduct(const std::vector<double>& vec_1, const std::vector<double>& vec_2)
{
  return InnerProduct(vec_1.size(), &vec_1[0], &vec_2[0]);
}

void
LU_Decompose(int num_entries, MatrixContainer& mat, int* index);

void
LU_Solve(int num_entries, MatrixContainer& mat, double* vec, int* index);

void
LU_SolveSystem(int num_entries, MatrixContainer& mat, double* vec, int* scratch);

struct CGScratchSpace
{
  CGScratchSpace() : len_(0) {}

  void
  Resize(int len)
  {
    if (len != len_) {
      len_ = len;
      d.resize(len_);
      r.resize(len_);
      s.resize(len_);
      q.resize(len_);
      M.AllocateDiagonalMatrix(len_);
    }
  }

  int                 len_;
  std::vector<double> d;
  std::vector<double> r;
  std::vector<double> s;
  std::vector<double> q;
  CRSMatrixContainer  M;
};

bool
CG_SolveSystem(
    nimble::CRSMatrixContainer& A,
    const double*               b,
    CGScratchSpace&             cg_scratch,
    double*                     x,
    int&                        num_iterations);

}  // namespace nimble

#endif
