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

#include <nimble_kokkos_block.h>
#include <nimble_kokkos_material_factory.h>

namespace nimble_kokkos {

void
Block::Initialize(std::string const& model_material_parameters, int num_elements, MaterialFactory& factory)
{
  model_material_parameters_ = model_material_parameters;
  InstantiateElement();
  int num_material_points = num_elements * element_->NumIntegrationPointsPerElement();
  InstantiateMaterialModel(num_material_points, factory);
}

void
Block::InstantiateMaterialModel(int num_material_points, MaterialFactory& factory)
{
  factory.parse_and_create(model_material_parameters_, num_material_points);
  material_        = factory.get_material_host();
  material_device_ = factory.get_material_device();
  ngp_lame_data_   = factory.get_ngp_lame_data();
}

void
Block::InstantiateElement()
{
  // instantiate the element on the host (eventually we won't need this)
  element_ = std::make_shared<nimble::HexElement>();

  // instantiate the element on the device
  element_device_ = static_cast<nimble::Element*>(Kokkos::kokkos_malloc<>("Element", sizeof(nimble::HexElement)));
  nimble::Element* pointer_that_lives_on_the_stack = element_device_;
  Kokkos::parallel_for(
      1, KOKKOS_LAMBDA(int) { new (pointer_that_lives_on_the_stack) nimble::HexElement(); });
}

}  // namespace nimble_kokkos
