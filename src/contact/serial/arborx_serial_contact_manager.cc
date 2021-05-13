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

#ifdef NIMBLE_HAVE_ARBORX

#include "arborx_serial_contact_manager.h"

#include "nimble_data_manager.h"
#include "nimble_defs.h"
#include "nimble_model_data_base.h"
#include "nimble_timer.h"
#include "nimble_vector_communicator.h"

#ifdef NIMBLE_HAVE_KOKKOS
#include "nimble_kokkos_model_data.h"
#endif

#include <ArborX.hpp>
#include <Kokkos_Core.hpp>
#include <iostream>
#include <random>
#include <vector>

#include "contact/arborx_utils.h"


namespace nimble {

using memory_space = nimble_kokkos::kokkos_device_memory_space;
using arborx_bvh   = ArborX::BVH<memory_space>;

/*!
 * Contact Manager specific to ArborX library
 *
 * @param interface Contact Interface
 */
ArborXSerialContactManager::ArborXSerialContactManager(
    std::shared_ptr<ContactInterface> interface,
    nimble::DataManager&              data_manager)
    : SerialContactManager{interface, data_manager}
{
  /// TODO Ask NM & RJ whether this is needed at constructor time
  Kokkos::View<int*, nimble_kokkos::kokkos_device> indices("indices", 0);
  Kokkos::View<int*, nimble_kokkos::kokkos_device> offset("offset", 0);
  updateCollisionData(indices, offset);
}

void
ArborXSerialContactManager::updateCollisionData(
    Kokkos::View<int*, nimble_kokkos::kokkos_device>& indices,
    Kokkos::View<int*, nimble_kokkos::kokkos_device>& offset)
{
  //
  // primitives are faces
  //
  arborx_bvh bvh{nimble_kokkos::kokkos_device_execution_space{}, contact_faces_d_};

  //
  // indices : position of the primitives that satisfy the predicates.
  // offsets : predicate offsets in indices.
  //
  // indices stores the indices of the objects that satisfy the predicates.
  // offset stores the locations in the indices view that start a predicate,
  // that is, predicates(q) is satisfied by indices(o) for primitives(q) <= o <
  // primitives(q+1).
  //
  // Following the usual convention, offset(n) = nnz, where n is the number of
  // queries that were performed and nnz is the total number of collisions.
  //
  // (From
  // https://github.com/arborx/ArborX/wiki/ArborX%3A%3ABoundingVolumeHierarchy%3A%3Aquery
  // )
  //

  // Number of queries, n = size of contact_nodes_d
  // Size of offset = n + 1
  bvh.query(nimble_kokkos::kokkos_device_execution_space{}, contact_nodes_d_, indices, offset);
}

template <typename ContactManagerType>
struct ArborXCallback {
  ContactManagerType contact_manager_;

  template <typename Query>
  KOKKOS_FUNCTION void operator()(Query const &query, int j) const {
    auto &myNode =
        contact_manager_.contact_nodes_d_(ArborX::getData(query));
    auto &myFace = contact_manager_.contact_faces_d_(j);

    double gap = 0.0;
    double normal[3] = {0., 0., 0.};
    bool inside = false;
    double facet_coordinates[3] = {0., 0., 0.};

    //--- Determine whether the node is projected inside the triangular face
    ContactManager::Projection(myNode, myFace, inside, gap, &normal[0],
                               &facet_coordinates[0]);
    if (inside) {
      myFace.set_contact_status(true);
      myNode.set_contact_status(true);
      contact_manager_.EnforceNodeFaceInteraction(
          myNode, myFace, gap, normal, facet_coordinates,
          contact_manager_.force_d_);
    }
  }
};

void
ArborXSerialContactManager::ComputeSerialContactForce(int step, bool debug_output, nimble::Viewify<2> contact_force)
{
  //--- Constraint per ContactManager::ComputeContactForce
  if (penalty_parameter_ <= 0.0) {
    throw std::logic_error("\nError in ComputeContactForce(), invalid penalty_parameter.\n");
  }

  if (model_data == nullptr) {
    auto model_ptr = data_manager_.GetMacroScaleData().get();
    model_data     = dynamic_cast<nimble_kokkos::ModelData*>(model_ptr);
  }

  auto field_ids      = data_manager_.GetFieldIDs();
  auto displacement_d = model_data->GetDeviceVectorNodeData(field_ids.displacement);

  auto contact_force_h = model_data->GetHostVectorNodeData(field_ids.contact_force);
  auto contact_force_d = model_data->GetDeviceVectorNodeData(field_ids.contact_force);
  Kokkos::deep_copy(contact_force_d, (double)(0.0));

  this->ApplyDisplacements(displacement_d);

  // Steps per (DJL)
  // entities are stored in contact_nodes_d, contact_nodes_h
  // 1) box box search
  // 2) node-face projection
  // 3) culling
  // 4) enforcement

  //--- Set vector to store force
  ContactManager::ZeroContactForce();

  // Reset the contact_status flags
  for (size_t jj = 0; jj < contact_faces_d_.extent(0); ++jj)
    contact_faces_d_(jj).set_contact_status(false);

  for (size_t jj = 0; jj < contact_nodes_d_.extent(0); ++jj)
    contact_nodes_d_(jj).set_contact_status(false);

  //--- Update the geometric collision information
  this->startTimer("Contact::EnforceInteraction");
  this->startTimer("ArborX::Search");
  arborx_bvh bvh{nimble_kokkos::kokkos_device_execution_space{},
                 contact_faces_d_};
  auto &contact_manager = *this;
  bvh.query(nimble_kokkos::kokkos_device_execution_space{}, contact_nodes_d_,
            ArborXCallback<decltype(contact_manager)>{contact_manager});
  this->stopTimer("ArborX::Search");
  this->stopTimer("Contact::EnforceInteraction");

  this->GetForces(contact_force_d);
  Kokkos::deep_copy(contact_force_h, contact_force_d);

}

}  // namespace nimble

#endif
