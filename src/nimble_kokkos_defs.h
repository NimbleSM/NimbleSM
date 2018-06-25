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

#ifndef NIMBLE_KOKKOS_DEFS_H
#define NIMBLE_KOKKOS_DEFS_H

#ifdef NIMBLE_HAVE_KOKKOS

#define NIMBLE_FUNCTION KOKKOS_FUNCTION
#define NIMBLE_INLINE_FUNCTION KOKKOS_INLINE_FUNCTION

#include "Kokkos_Core.hpp"

// Define Memory Space
#ifdef KOKKOS_ENABLE_CUDA
  using kokkos_device_memory_space = Kokkos::CudaSpace; // Kokkos::HostSpace, Kokkos::CudaSpace, Kokkos::CudaUVMSpace, Kokkos::CudaHostPinnedSpace, etc.
#else
  using kokkos_device_memory_space = Kokkos::HostSpace;
#endif
using kokkos_host_memory_space = Kokkos::HostSpace;

namespace nimble_kokkos {

// Layout should not be specified unless you're using teams (and maybe not even if you are)
using Layout = kokkos_device_memory_space::execution_space::array_layout;

enum class FieldType : int
{
  HostScalar, DeviceScalar
, HostVector, DeviceVector
, HostSymTensor, DeviceSymTensor
, HostFullTensor, DeviceFullTensor
};

inline
std::ostream & operator << ( std::ostream & out, const FieldType f )
{
  switch( f )
  {
  case FieldType::HostScalar:
    out << "HostScalar";
    break;
  case FieldType::DeviceScalar:
    out << "DeviceScalar";
    break;
  case FieldType::HostVector:
    out << "HostVector";
    break;
  case FieldType::DeviceVector:
    out << "DeviceVector";
    break;
  case FieldType::HostSymTensor:
    out << "HostSymTensor";
    break;
  case FieldType::DeviceSymTensor:
    out << "DeviceSymTensor";
    break;
  case FieldType::HostFullTensor:
    out << "HostFullTensor";
    break;
  case FieldType::DeviceFullTensor:
    out << "DeviceFullTensor";
    break;
  };
  return out;
}

class FieldBase
{
public:
  FieldType type() const noexcept { return type_; }

  const std::string& name() const noexcept { return name_; }

  FieldBase( const std::string & arg_name
           , const FieldType arg_type
           )
    : type_{arg_type}
    , name_{arg_name}
  {}

  FieldBase() = default;

  FieldBase( const FieldBase &  ) = default;
  FieldBase(       FieldBase && ) = default;

  FieldBase & operator=( const FieldBase &  ) = default;
  FieldBase & operator=(       FieldBase && ) = default;

  virtual ~FieldBase() = default;

  static constexpr int NUM_NODES_IN_HEX = 8;
  static constexpr int NUM_INTEGRATION_POINTS_IN_HEX = 8;

private:
  FieldType type_   {};
  std::string name_ {};
};

template < FieldType FType >
class Field;

template <>
class Field< FieldType::HostScalar >
  : public FieldBase
{
public:

  using View = Kokkos::View< double *, kokkos_host_memory_space::execution_space >;

  Field( const std::string & name, int num_entries )
    : FieldBase( name, FieldType::HostScalar )
    , data_( name, num_entries )
  {}

  Field( View const & v )
    : FieldBase( v.label(), FieldType::HostScalar )
    , data_(v)
  {}

  View data() const noexcept { return data_; }

private:
  View data_;

};

template <>
class Field< FieldType::DeviceScalar >
  : public FieldBase
{
public:

  using View = Kokkos::View< double *, Layout, kokkos_device_memory_space::execution_space >;
  using AtomicView = Kokkos::View< double *, Layout, kokkos_device_memory_space::execution_space, Kokkos::MemoryTraits<Kokkos::Atomic> >;
  using GatheredView = Kokkos::View< double *[NUM_NODES_IN_HEX], Layout, kokkos_device_memory_space::execution_space >;
  using GatheredSubView = decltype(Kokkos::subview(*(GatheredView*)(0), (int)(0), Kokkos::ALL));

 Field( const std::string & name, int num_entries )
    : FieldBase( name, FieldType::DeviceScalar )
    , data_( name, num_entries )
  {}

  Field( View const & v )
    : FieldBase( v.label(), FieldType::DeviceScalar )
    , data_(v)
  {}

  View data() const noexcept { return data_; }

private:
  View data_;

};

template <>
class Field< FieldType::HostVector >
  : public FieldBase
{
public:

  using View = Kokkos::View< double *[3], Layout, kokkos_host_memory_space::execution_space >;

  Field( const std::string & name
       , int num_entries
       )
    : FieldBase( name, FieldType::HostVector )
    , data_( name, num_entries )
  {}

  Field( View const & v )
    : FieldBase( v.label(), FieldType::HostVector )
    , data_(v)
  {}

  View data() const noexcept { return data_; }

private:
  View data_;
};

template <>
class Field< FieldType::DeviceVector >
  : public FieldBase
{
public:

  using View = Kokkos::View< double *[3], Layout, kokkos_device_memory_space::execution_space >;
  using AtomicView = Kokkos::View< double *[3], Layout, kokkos_device_memory_space::execution_space, Kokkos::MemoryTraits<Kokkos::Atomic> >;
  using GatheredView = Kokkos::View< double *[NUM_NODES_IN_HEX][3], Layout, kokkos_device_memory_space::execution_space >;
  using GatheredSubView = decltype(Kokkos::subview(*(GatheredView*)(0), (int)(0), Kokkos::ALL, Kokkos::ALL));

  Field( const std::string & name
       , int num_entries
       )
    : FieldBase( name, FieldType::DeviceVector )
    , data_( name, num_entries )
  {}

  Field( View const & v )
    : FieldBase( v.label(), FieldType::DeviceVector )
    , data_(v)
  {}

  View data() const noexcept { return data_; }

private:
  View data_;
};


template <>
class Field< FieldType::HostSymTensor >
  : public FieldBase
{
public:

  using View = Kokkos::View< double *[NUM_INTEGRATION_POINTS_IN_HEX][6], Layout, kokkos_host_memory_space::execution_space >;

  Field( const std::string & name
       , int num_entries
       )
    : FieldBase( name, FieldType::HostSymTensor )
    , data_( name, num_entries )
  {}

  Field( View const & v )
    : FieldBase( v.label(), FieldType::HostSymTensor )
    , data_(v)
  {}


  View data() const noexcept { return data_; }

private:
  View data_;
};

template <>
class Field< FieldType::DeviceSymTensor >
  : public FieldBase
{
public:

  using View = Kokkos::View< double *[NUM_INTEGRATION_POINTS_IN_HEX][6], Layout, kokkos_device_memory_space::execution_space >;
  using SubView = decltype(Kokkos::subview(*(View*)(0), (int)(0), Kokkos::ALL, Kokkos::ALL));
  using SingleEntryView = decltype(Kokkos::subview(*(View*)(0), (int)(0), (int)(0), Kokkos::ALL));

 Field( const std::string & name
       , int num_entries
       )
    : FieldBase( name, FieldType::DeviceSymTensor )
    , data_( name, num_entries )
  {}

  Field( View const & v )
    : FieldBase( v.label(), FieldType::DeviceSymTensor )
    , data_(v)
  {}

  View data() const noexcept { return data_; }

private:
  View data_;
};


template <>
class Field< FieldType::HostFullTensor >
  : public FieldBase
{
public:

  using View = Kokkos::View< double *[NUM_INTEGRATION_POINTS_IN_HEX][9], Layout, kokkos_host_memory_space::execution_space >;

  Field( const std::string & name
       , int num_entries
       )
    : FieldBase( name, FieldType::HostFullTensor )
    , data_( name, num_entries )
  {}

  Field( View const & v )
    : FieldBase( v.label(), FieldType::HostFullTensor )
    , data_(v)
  {}

  View data() const noexcept { return data_; }

private:
  View data_;
};

template <>
class Field< FieldType::DeviceFullTensor >
  : public FieldBase
{
public:

  using View = Kokkos::View< double *[NUM_INTEGRATION_POINTS_IN_HEX][9], Layout, kokkos_device_memory_space::execution_space >;
  using SubView = decltype(Kokkos::subview(*(View*)(0), (int)(0), Kokkos::ALL, Kokkos::ALL));
  using SingleEntryView = decltype(Kokkos::subview(*(View*)(0), (int)(0), (int)(0), Kokkos::ALL));

  Field( const std::string & name
       , int num_entries
       )
    : FieldBase( name, FieldType::DeviceFullTensor )
    , data_( name, num_entries )
  {}

  Field( View const & v )
    : FieldBase( v.label(), FieldType::DeviceFullTensor )
    , data_(v)
  {}

  View data() const noexcept { return data_; }

private:
  View data_;
};

// Field< FielType::HostVector>::View is a Kokkos::View< double *[3], Layout, kokkos_host_memory_space::execution_space >
typedef Field< FieldType::HostScalar >::View                                          HostScalarView;
typedef Field< FieldType::HostVector >::View                                          HostVectorView;
typedef Kokkos::View< int *, Layout, kokkos_host_memory_space::execution_space >      HostElementConnectivityView; // TODO THIS SHOULD BE A 2D ARRAY, BUT IT'S TRICKY BECAUSE NUM NODES PER ELEMENT IS NOT KNOWN

typedef Field< FieldType::DeviceScalar >::View                                        DeviceScalarView;
typedef Field< FieldType::DeviceScalar >::GatheredView                                DeviceScalarGatheredView;
typedef Field< FieldType::DeviceScalar >::GatheredSubView                             DeviceScalarGatheredSubView;
typedef Field< FieldType::DeviceVector >::View                                        DeviceVectorView;
typedef Field< FieldType::DeviceVector >::GatheredView                                DeviceVectorGatheredView;
typedef Field< FieldType::DeviceVector >::GatheredSubView                             DeviceVectorGatheredSubView;
typedef Field< FieldType::DeviceFullTensor >::View                                    DeviceFullTensorView;
typedef Field< FieldType::DeviceFullTensor >::SubView                                 DeviceFullTensorSubView;
typedef Field< FieldType::DeviceFullTensor >::SingleEntryView                         DeviceFullTensorSingleEntryView;
typedef Field< FieldType::DeviceSymTensor >::View                                     DeviceSymTensorView;
typedef Field< FieldType::DeviceSymTensor >::SubView                                  DeviceSymTensorSubView;
typedef Field< FieldType::DeviceSymTensor >::SingleEntryView                          DeviceSymTensorSingleEntryView;
typedef Kokkos::View< int *, Layout, kokkos_device_memory_space::execution_space >    DeviceElementConnectivityView;
}

// else for NIMBLE_HAVE_KOKKOS
#else

#define NIMBLE_FUNCTION
#define NIMBLE_INLINE_FUNCTION inline

// endif for NIMBLE_HAVE_KOKKOS
#endif

// endif for NIMBLE_KOKKOS_DEFS_H
#endif
