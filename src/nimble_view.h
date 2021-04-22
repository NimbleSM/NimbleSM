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

#include <array>
#include <memory>
#include <type_traits>

#include "nimble_defs.h"

enum class FieldEnum : int
{
  Scalar     = 1,
  Vector     = 3,
  SymTensor  = 6,
  FullTensor = 9
};

namespace nimble {

namespace details {

template <std::size_t N>
class AXPYResult;

}

template <std::size_t N = 2, class Scalar = double>
class Viewify
{
 public:
  Viewify() : data_(nullptr)
  {
    len_.fill(0);
    stride_.fill(0);
  }

  Viewify(Scalar* data, std::array<int, N> len, std::array<int, N> stride)
      : data_(data), len_(len), stride_(stride)
  {
  }

  template <std::size_t NN = N, typename = typename std::enable_if<(NN == 1)>::type >
  Viewify(Scalar* data, int len)
      : data_(data), len_({len}), stride_({1})
  { }

  template <std::size_t NN = N>
  NIMBLE_INLINE_FUNCTION typename std::enable_if<(NN == 1), Scalar>::type&
  operator()(int i)
  {
    return data_[i];
  }

  template <std::size_t NN = N>
  NIMBLE_INLINE_FUNCTION typename std::enable_if<(NN == 1), const Scalar>::type
  operator()(int i) const
  {
    return data_[i];
  }

  template <std::size_t NN = N>
  NIMBLE_INLINE_FUNCTION typename std::enable_if<(NN == 2), Scalar>::type&
  operator()(int i, int j)
  {
    return data_[i * stride_[0] + j * stride_[1]];
  }

  template <std::size_t NN = N>
  NIMBLE_INLINE_FUNCTION typename std::enable_if<(NN == 2), const Scalar>::type
  operator()(int i, int j) const
  {
    return data_[i * stride_[0] + j * stride_[1]];
  }

  void
  zero()
  {
    int mySize = stride_[0] * len_[0];
    for (int ii = 0; ii < mySize; ++ii) data_[ii] = (Scalar)(0);
  }

  void
  copy(const Viewify<N>& ref)
  {
    int mySize = stride_[0] * len_[0];
    for (int ii = 0; ii < mySize; ++ii) data_[ii] = ref.data_[ii];
  }

  Scalar*
  data() const
  {
    return data_;
  }

  template< class T = Scalar, typename = typename std::enable_if< !std::is_const<T>::value >::type >
  Viewify<N, Scalar>&
  operator+=(const details::AXPYResult<N>& rhs);

  std::array<int, N>
  size() const
  {
    return len_;
  }

  std::array<int, N>
  stride() const
  {
    return stride_;
  }

 protected:
  Scalar*            data_;
  std::array<int, N> len_;
  std::array<int, N> stride_;
};

//----------------------------------

namespace details {

template <std::size_t N>
struct AXPYResult
{
  AXPYResult(const nimble::Viewify<N> A, double alpha) : A_(A), alpha_(alpha) {}

  /// \brief Compute "dest = (destCoef) * dest + rhsCoef * alpha_ * A_"
  ///
  /// \param[in] destCoef Scalar
  /// \param[in] rhsCoef Scalar
  /// \param[out] dest Viewify "vector" storing the result
  ///
  /// \note The routine does not check whether A_ and dest have the same length
  /// and the same stride.
  void
  assignTo(double destCoef, double rhsCoef, nimble::Viewify<N>& dest) const;

  const nimble::Viewify<N> A_;
  double                   alpha_;
};

template <std::size_t N>
void
AXPYResult<N>::assignTo(double destCoef, double rhsCoef, nimble::Viewify<N>& dest) const
{
  const double prod = alpha_ * rhsCoef;

  auto       len    = dest.size();
  auto       stride = dest.stride();
  const long isize  = len[0] * stride[0];

  double* data     = dest.data();
  double* rhs_data = A_.data();

  if (destCoef == 1.0) {
    for (long ii = 0; ii < isize; ++ii) data[ii] += prod * rhs_data[ii];
  } else if (destCoef == 0.0) {
    for (long ii = 0; ii < isize; ++ii) data[ii] = prod * rhs_data[ii];
  } else {
    for (long ii = 0; ii < isize; ++ii) data[ii] = destCoef * data[ii] + prod * rhs_data[ii];
  }
}

}  // namespace details

//----------------------------------

template <std::size_t N, class Scalar>
template <class T, typename >
Viewify<N, Scalar>&
Viewify<N, Scalar>::operator+=(const details::AXPYResult<N>& rhs)
{
  rhs.assignTo(1.0, 1.0, *this);
  return *this;
}

template <std::size_t N>
details::AXPYResult<N>
operator*(double alpha, const nimble::Viewify<N>& A)
{
  return {A, alpha};
}

//----------------------------------

template <FieldEnum FieldT>
class View
{
 public:
  explicit View(int num_entries) : num_entries_(num_entries)
  {
    data_ = std::shared_ptr<double>(new double[num_entries_ * static_cast<int>(FieldT)], [](double* p) { delete[] p; });
    data_ptr_ = data_.get();
  }

  ~View()               = default;
  View(const View&)     = default;
  View(View&&) noexcept = default;
  View&
  operator=(const View&) = default;
  View&
  operator=(View&&) noexcept = default;

  double&
  operator()(int i_entry) const
  {
    static_assert(FieldT == FieldEnum::Scalar, "Operator(int i_entry) called for non-scalar data.");
    return data_ptr_[i_entry];
  }

  double&
  operator()(int i_entry, int i_coord) const
  {
    static_assert(FieldT != FieldEnum::Scalar, "Operator(int i_entry, int i_coord) called for scalar data.");
    return data_ptr_[i_entry * static_cast<int>(FieldT) + i_coord];
  }

  int
  rank()
  {
    if (FieldT == FieldEnum::Scalar) { return 1; }
    return 2;
  }

  int
  extent(int dim)
  {
    if (dim == 0) {
      return num_entries_;
    } else if (dim == 1) {
      return static_cast<int>(FieldT);
    } else {
      return 1;
    }
  }

#ifdef NIMBLE_HAVE_DARMA
  void
  serialize(darma_runtime::serialization::SimplePackUnpackArchive& ar)
  {
    // The purpose for the serialize call can be determined with:
    // ar.is_sizing()
    // ar.is_packing()
    // ar.is_unpacking()

    const int           size = num_entries_ * static_cast<int>(FieldT);
    std::vector<double> data_vec(size);

    if (!ar.is_unpacking()) {
      for (int i = 0; i < size; i++) { data_vec[i] = data_ptr_[i]; }
    }

    ar | num_entries_ | data_vec;

    if (ar.is_unpacking()) {
      data_ =
          std::shared_ptr<double>(new double[num_entries_ * static_cast<int>(FieldT)], [](double* p) { delete[] p; });
      data_ptr_ = data_.get();
      for (int i = 0; i < size; i++) { data_ptr_[i] = data_vec[i]; }
    }
  }
#endif

 private:
  int                     num_entries_;
  std::shared_ptr<double> data_;
  double*                 data_ptr_;
};

}  // namespace nimble

#endif
