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

namespace nimble_kokkos {

// Define HOST execution space, memory space, and device
using kokkos_host_execution_space = Kokkos::Serial::execution_space;
#ifdef KOKKOS_ENABLE_CUDA_UVM
  using kokkos_host_mirror_memory_space = Kokkos::CudaUVMSpace::memory_space;
#else
  using kokkos_host_mirror_memory_space = Kokkos::Serial::memory_space;
#endif
using kokkos_host = Kokkos::Device<kokkos_host_execution_space, kokkos_host_mirror_memory_space>;

// Define DEVICE execution space, memory space, and device
#ifdef KOKKOS_ENABLE_CUDA_UVM
  using kokkos_device_execution_space = Kokkos::Cuda::execution_space;
  using kokkos_device_memory_space = Kokkos::CudaUVMSpace::memory_space;
#else
  #ifdef KOKKOS_ENABLE_CUDA
    using kokkos_device_execution_space = Kokkos::Cuda::execution_space;
    using kokkos_device_memory_space = Kokkos::Cuda::memory_space;
  #else
    using kokkos_device_execution_space = Kokkos::Serial::execution_space;
    using kokkos_device_memory_space = Kokkos::Serial::memory_space;
  #endif
#endif
using kokkos_device = Kokkos::Device<kokkos_device_execution_space, kokkos_device_memory_space>;

using kokkos_layout = kokkos_device_execution_space::array_layout;

// Layout should not be specified unless you're using teams (and maybe not even if you are)
//using Layout = kokkos_device_memory_space::execution_space::array_layout;

enum class FieldType : int
{
  HostScalarNode, DeviceScalarNode
, HostVectorNode, DeviceVectorNode
, HostSymTensorIntPt, DeviceSymTensorIntPt
, HostFullTensorIntPt, DeviceFullTensorIntPt
, HostScalarElem, DeviceScalarElem
, HostSymTensorElem, DeviceSymTensorElem
, HostFullTensorElem, DeviceFullTensorElem
};

inline
std::ostream & operator << ( std::ostream & out, const FieldType f )
{
  switch( f )
  {
  case FieldType::HostScalarNode:
    out << "HostScalarNode";
    break;
  case FieldType::DeviceScalarNode:
    out << "DeviceScalarNode";
    break;
  case FieldType::HostVectorNode:
    out << "HostVectorNode";
    break;
  case FieldType::DeviceVectorNode:
    out << "DeviceVectorNode";
    break;
  case FieldType::HostSymTensorIntPt:
    out << "HostSymTensorIntPt";
    break;
  case FieldType::DeviceSymTensorIntPt:
    out << "DeviceSymTensorIntPt";
    break;
  case FieldType::HostFullTensorIntPt:
    out << "HostFullTensorIntPt";
    break;
  case FieldType::DeviceFullTensorIntPt:
    out << "DeviceFullTensorIntPt";
    break;
  case FieldType::HostScalarElem:
    out << "HostScalarElem";
    break;
  case FieldType::DeviceScalarElem:
    out << "DeviceScalarElem";
    break;
  case FieldType::HostSymTensorElem:
    out << "HostSymTensorElem";
    break;
  case FieldType::DeviceSymTensorElem:
    out << "DeviceSymTensorElem";
    break;
  case FieldType::HostFullTensorElem:
    out << "HostFullTensorElem";
    break;
  case FieldType::DeviceFullTensorElem:
    out << "DeviceFullTensorElem";
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
class Field< FieldType::HostScalarNode >
  : public FieldBase
{
public:

  // (node)
  using View = Kokkos::View< double *, kokkos_host >;

  Field( const std::string & name, int num_entries )
    : FieldBase( name, FieldType::HostScalarNode )
    , data_( name, num_entries )
  {}

  Field( View const & v )
    : FieldBase( v.label(), FieldType::HostScalarNode )
    , data_(v)
  {}

  View data() const noexcept { return data_; }

private:
  View data_;

};

template <>
class Field< FieldType::DeviceScalarNode >
  : public FieldBase
{
public:

  // (node)
  using View = Kokkos::View< double *, kokkos_layout, kokkos_device >;
  using AtomicView = Kokkos::View< double *, kokkos_layout, kokkos_device, Kokkos::MemoryTraits<Kokkos::Atomic> >;
  using GatheredView = Kokkos::View< double *[NUM_NODES_IN_HEX], kokkos_layout, kokkos_device >;
  using GatheredSubView = decltype(Kokkos::subview(*(GatheredView*)(0), (int)(0), Kokkos::ALL));

 Field( const std::string & name, int num_entries )
    : FieldBase( name, FieldType::DeviceScalarNode )
    , data_( name, num_entries )
  {}

  Field( View const & v )
    : FieldBase( v.label(), FieldType::DeviceScalarNode )
    , data_(v)
  {}

  View data() const noexcept { return data_; }

private:
  View data_;

};

template <>
class Field< FieldType::HostVectorNode >
  : public FieldBase
{
public:

  // (node, coordinate)
  using View = Kokkos::View< double *[3], kokkos_layout, kokkos_host >;

  Field( const std::string & name
       , int num_entries
       )
    : FieldBase( name, FieldType::HostVectorNode )
    , data_( name, num_entries )
  {}

  Field( View const & v )
    : FieldBase( v.label(), FieldType::HostVectorNode )
    , data_(v)
  {}

  View data() const noexcept { return data_; }

private:
  View data_;
};

template <>
class Field< FieldType::DeviceVectorNode >
  : public FieldBase
{
public:

  // (node, coordinate)
  using View = Kokkos::View< double *[3], kokkos_layout, kokkos_device >;
  using AtomicView = Kokkos::View< double *[3], kokkos_layout, kokkos_device, Kokkos::MemoryTraits<Kokkos::Atomic> >;
  using GatheredView = Kokkos::View< double *[NUM_NODES_IN_HEX][3], kokkos_layout, kokkos_device >;
  using GatheredSubView = decltype(Kokkos::subview(*(GatheredView*)(0), (int)(0), Kokkos::ALL, Kokkos::ALL));

  Field( const std::string & name
       , int num_entries
       )
    : FieldBase( name, FieldType::DeviceVectorNode )
    , data_( name, num_entries )
  {}

  Field( View const & v )
    : FieldBase( v.label(), FieldType::DeviceVectorNode )
    , data_(v)
  {}

  View data() const noexcept { return data_; }

private:
  View data_;
};

template <>
class Field< FieldType::HostSymTensorIntPt >
  : public FieldBase
{
public:

  // (elem, ipt, tensor_index)
  using View = Kokkos::View< double *[NUM_INTEGRATION_POINTS_IN_HEX][6], kokkos_layout, kokkos_host >;

  Field( const std::string & name
       , int num_entries
       )
    : FieldBase( name, FieldType::HostSymTensorIntPt )
    , data_( name, num_entries )
  {}

  Field( View const & v )
    : FieldBase( v.label(), FieldType::HostSymTensorIntPt )
    , data_(v)
  {}


  View data() const noexcept { return data_; }

private:
  View data_;
};

template <>
class Field< FieldType::DeviceSymTensorIntPt >
  : public FieldBase
{
public:

  // (elem, ipt, tensor_index)
  using View = Kokkos::View< double *[NUM_INTEGRATION_POINTS_IN_HEX][6], kokkos_layout, kokkos_device >;
  using SubView = decltype(Kokkos::subview(*(View*)(0), (int)(0), Kokkos::ALL, Kokkos::ALL));
  using SingleEntryView = decltype(Kokkos::subview(*(View*)(0), (int)(0), (int)(0), Kokkos::ALL));

 Field( const std::string & name
       , int num_entries
       )
    : FieldBase( name, FieldType::DeviceSymTensorIntPt )
    , data_( name, num_entries )
  {}

  Field( View const & v )
    : FieldBase( v.label(), FieldType::DeviceSymTensorIntPt )
    , data_(v)
  {}

  View data() const noexcept { return data_; }

private:
  View data_;
};

template <>
class Field< FieldType::HostFullTensorIntPt >
  : public FieldBase
{
public:

  // (elem, ipt, tensor_index)
  using View = Kokkos::View< double *[NUM_INTEGRATION_POINTS_IN_HEX][9], kokkos_layout, kokkos_host >;

  Field( const std::string & name
       , int num_entries
       )
    : FieldBase( name, FieldType::HostFullTensorIntPt )
    , data_( name, num_entries )
  {}

  Field( View const & v )
    : FieldBase( v.label(), FieldType::HostFullTensorIntPt )
    , data_(v)
  {}

  View data() const noexcept { return data_; }

private:
  View data_;
};

template <>
class Field< FieldType::DeviceFullTensorIntPt >
  : public FieldBase
{
public:

  // (elem, ipt, tensor_index)
  using View = Kokkos::View< double *[NUM_INTEGRATION_POINTS_IN_HEX][9], kokkos_layout, kokkos_device >;
  using SubView = decltype(Kokkos::subview(*(View*)(0), (int)(0), Kokkos::ALL, Kokkos::ALL));
  using SingleEntryView = decltype(Kokkos::subview(*(View*)(0), (int)(0), (int)(0), Kokkos::ALL));

  Field( const std::string & name
       , int num_entries
       )
    : FieldBase( name, FieldType::DeviceFullTensorIntPt )
    , data_( name, num_entries )
  {}

  Field( View const & v )
    : FieldBase( v.label(), FieldType::DeviceFullTensorIntPt )
    , data_(v)
  {}

  View data() const noexcept { return data_; }

private:
  View data_;
};

template <>
class Field< FieldType::HostScalarElem >
  : public FieldBase
{
public:

  // (elem)
  using View = Kokkos::View< double *, kokkos_host >;

  Field( const std::string & name, int num_entries )
    : FieldBase( name, FieldType::HostScalarElem )
    , data_( name, num_entries )
  {}

  Field( View const & v )
    : FieldBase( v.label(), FieldType::HostScalarElem )
    , data_(v)
  {}

  View data() const noexcept { return data_; }

private:
  View data_;

};

template <>
class Field< FieldType::DeviceScalarElem >
  : public FieldBase
{
public:

  // (node)
  using View = Kokkos::View< double *, kokkos_layout, kokkos_device >;
  using SingleEntryView = decltype(Kokkos::subview(*(View*)(0), (int)(0)));

 Field( const std::string & name, int num_entries )
    : FieldBase( name, FieldType::DeviceScalarElem )
    , data_( name, num_entries )
  {}

  Field( View const & v )
    : FieldBase( v.label(), FieldType::DeviceScalarElem )
    , data_(v)
  {}

  View data() const noexcept { return data_; }

private:
  View data_;

};

template <>
class Field< FieldType::HostFullTensorElem >
  : public FieldBase
{
public:

  // (elem, ipt, tensor_index)
  using View = Kokkos::View< double *[9], kokkos_layout, kokkos_host >;

  Field( const std::string & name
       , int num_entries
       )
    : FieldBase( name, FieldType::HostFullTensorElem )
    , data_( name, num_entries )
  {}

  Field( View const & v )
    : FieldBase( v.label(), FieldType::HostFullTensorElem )
    , data_(v)
  {}

  View data() const noexcept { return data_; }

private:
  View data_;
};

template <>
class Field< FieldType::DeviceFullTensorElem >
  : public FieldBase
{
public:

  // (elem, ipt, tensor_index)
  using View = Kokkos::View< double *[9], kokkos_layout ,kokkos_device >;
  using SingleEntryView = decltype(Kokkos::subview(*(View*)(0), (int)(0), Kokkos::ALL));

  Field( const std::string & name
       , int num_entries
       )
    : FieldBase( name, FieldType::DeviceFullTensorElem )
    , data_( name, num_entries )
  {}

  Field( View const & v )
    : FieldBase( v.label(), FieldType::DeviceFullTensorElem )
    , data_(v)
  {}

  View data() const noexcept { return data_; }

private:
  View data_;
};

template <>
class Field< FieldType::HostSymTensorElem >
  : public FieldBase
{
public:

  // (elem, ipt, tensor_index)
  using View = Kokkos::View< double *[6], kokkos_layout, kokkos_host >;

  Field( const std::string & name
       , int num_entries
       )
    : FieldBase( name, FieldType::HostSymTensorElem )
    , data_( name, num_entries )
  {}

  Field( View const & v )
    : FieldBase( v.label(), FieldType::HostSymTensorElem )
    , data_(v)
  {}


  View data() const noexcept { return data_; }

private:
  View data_;
};

template <>
class Field< FieldType::DeviceSymTensorElem >
  : public FieldBase
{
public:

  // (elem, ipt, tensor_index)
  using View = Kokkos::View< double *[6], kokkos_layout, kokkos_device >;
  using SingleEntryView = decltype(Kokkos::subview(*(View*)(0), (int)(0), Kokkos::ALL));

 Field( const std::string & name
       , int num_entries
       )
    : FieldBase( name, FieldType::DeviceSymTensorElem )
    , data_( name, num_entries )
  {}

  Field( View const & v )
    : FieldBase( v.label(), FieldType::DeviceSymTensorElem )
    , data_(v)
  {}

  View data() const noexcept { return data_; }

private:
  View data_;
};

typedef Field< FieldType::HostScalarNode >::View                                      HostScalarNodeView;
typedef Field< FieldType::HostVectorNode >::View                                      HostVectorNodeView;
typedef Field< FieldType::HostFullTensorIntPt >::View                                 HostFullTensorIntPtView;
typedef Field< FieldType::HostSymTensorIntPt >::View                                  HostSymTensorIntPtView;
typedef Field< FieldType::HostScalarElem >::View                                      HostScalarElemView;
typedef Field< FieldType::HostFullTensorElem >::View                                  HostFullTensorElemView;
typedef Field< FieldType::HostSymTensorElem >::View                                   HostSymTensorElemView;
typedef Kokkos::View< int*, kokkos_host >                                             HostIntegerArrayView;
typedef Kokkos::View< int*, kokkos_host >                                             HostElementConnectivityView; // TODO THIS SHOULD BE A 2D ARRAY, BUT IT'S TRICKY BECAUSE NUM NODES PER ELEMENT IS NOT KNOWN

typedef Field< FieldType::DeviceScalarNode >::View                                    DeviceScalarNodeView;
typedef Field< FieldType::DeviceScalarNode >::GatheredView                            DeviceScalarNodeGatheredView;
typedef Field< FieldType::DeviceScalarNode >::GatheredSubView                         DeviceScalarNodeGatheredSubView;
typedef Field< FieldType::DeviceVectorNode >::View                                    DeviceVectorNodeView;
typedef Field< FieldType::DeviceVectorNode >::GatheredView                            DeviceVectorNodeGatheredView;
typedef Field< FieldType::DeviceVectorNode >::GatheredSubView                         DeviceVectorNodeGatheredSubView;
typedef Field< FieldType::DeviceFullTensorIntPt >::View                               DeviceFullTensorIntPtView;
typedef Field< FieldType::DeviceFullTensorIntPt >::SubView                            DeviceFullTensorIntPtSubView;
typedef Field< FieldType::DeviceFullTensorIntPt >::SingleEntryView                    DeviceFullTensorIntPtSingleEntryView;
typedef Field< FieldType::DeviceSymTensorIntPt >::View                                DeviceSymTensorIntPtView;
typedef Field< FieldType::DeviceSymTensorIntPt >::SubView                             DeviceSymTensorIntPtSubView;
typedef Field< FieldType::DeviceSymTensorIntPt >::SingleEntryView                     DeviceSymTensorIntPtSingleEntryView;
typedef Field< FieldType::DeviceScalarElem >::View                                    DeviceScalarElemView;
typedef Field< FieldType::DeviceScalarElem >::SingleEntryView                         DeviceScalarElemSingleEntryView;
typedef Field< FieldType::DeviceFullTensorElem >::View                                DeviceFullTensorElemView;
typedef Field< FieldType::DeviceFullTensorElem >::SingleEntryView                     DeviceFullTensorElemSingleEntryView;
typedef Field< FieldType::DeviceSymTensorElem >::View                                 DeviceSymTensorElemView;
typedef Field< FieldType::DeviceSymTensorElem >::SingleEntryView                      DeviceSymTensorElemSingleEntryView;
typedef Kokkos::View< int*, kokkos_device >                                           DeviceIntegerArrayView;
typedef Kokkos::View< int*, kokkos_device >                                           DeviceElementConnectivityView;
}

// else for NIMBLE_HAVE_KOKKOS
#else

#define NIMBLE_FUNCTION
#define NIMBLE_INLINE_FUNCTION inline

// endif for NIMBLE_HAVE_KOKKOS
#endif

// endif for NIMBLE_KOKKOS_DEFS_H
#endif
