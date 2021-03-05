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

#ifndef NIMBLE_VIEW_H
#define NIMBLE_VIEW_H

#include <memory>

#include "nimble_defs.h"

class Viewify {

 public:

  Viewify() : data_(nullptr), dim_(0)
  {} // for NIMBLE_HAVE_UQ?

  Viewify(double * const data, int dim)
      : data_(data), dim_(dim) {}

  NIMBLE_INLINE_FUNCTION
  double& operator()(int i, int j) {
    return data_[i*dim_ + j];
  }

  NIMBLE_INLINE_FUNCTION
  const double& operator()(int i, int j) const {
    return data_[i*dim_ + j];
  }

private:

  double* const data_;
  int dim_;
};

enum class FieldEnum : int
{
  Scalar = 1
, Vector = 3
, SymTensor = 6
, FullTensor = 9
};

namespace nimble {

  template <FieldEnum FieldT>
  class View {

  public:

    View(int num_entries) : num_entries_(num_entries) {
      data_ = std::shared_ptr<double>(new double[num_entries_ * static_cast<int>(FieldT)], [](double* p){ delete[] p; });
      data_ptr_ = data_.get();
    }

    ~View() = default;
    View(const View &) = default;
    View(View &&) = default;
    View & operator=(const View &) = default;
    View & operator=(View &&) = default;

    double & operator() (int i_entry) const {
      static_assert(FieldT == FieldEnum::Scalar, "Operator(int i_entry) called for non-scalar data.");
      return data_ptr_[i_entry];
    }

    double & operator() (int i_entry, int i_coord) const {
      static_assert(FieldT != FieldEnum::Scalar, "Operator(int i_entry, int i_coord) called for scalar data.");
      return data_ptr_[i_entry * static_cast<int>(FieldT) + i_coord];
    }

    int rank() {
      if (FieldT == FieldEnum::Scalar) {
        return 1;
      }
      return 2;
    }

    int extent(int dim) {
      if (dim == 0) {
        return num_entries_;
      }
      else if (dim == 1) {
        return static_cast<int>(FieldT);
      }
      else {
        return 1;
      }
    }

#ifdef NIMBLE_HAVE_DARMA
  void serialize(darma_runtime::serialization::SimplePackUnpackArchive& ar) {

    // The purpose for the serialize call can be determined with:
    // ar.is_sizing()
    // ar.is_packing()
    // ar.is_unpacking()

    const int size = num_entries_ * static_cast<int>(FieldT);
    std::vector<double> data_vec(size);

    if (!ar.is_unpacking()) {
      for (int i=0 ; i<size ; i++) {
        data_vec[i] = data_ptr_[i];
      }
    }

    ar | num_entries_ | data_vec ;

    if (ar.is_unpacking()) {
      data_ = std::shared_ptr<double>(new double[num_entries_ * static_cast<int>(FieldT)], [](double* p){ delete[] p; });
      data_ptr_ = data_.get();
      for (int i=0 ; i<size ; i++) {
        data_ptr_[i] = data_vec[i];
      }
    }
  }
#endif

  private:

    int num_entries_;
    std::shared_ptr<double> data_;
    double * data_ptr_;
  };

}

#endif
