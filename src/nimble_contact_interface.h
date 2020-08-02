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

#ifndef SRC_NIMBLE_CONTACT_INTERFACE_H_
#define SRC_NIMBLE_CONTACT_INTERFACE_H_

#include <iostream>

#ifdef NIMBLE_HAVE_KOKKOS
  #include <nimble_kokkos_contact_defs.h>
  #include <nimble_kokkos_defs.h>
#endif

namespace nimble {

#ifdef NIMBLE_HAVE_KOKKOS

struct PenaltyContactEnforcement {
  PenaltyContactEnforcement() : penalty(0.0) {}

  KOKKOS_FORCEINLINE_FUNCTION
  void EnforceContact(ContactEntity &node, ContactEntity &face, int numNodeFaces, const double gap,
                      const double direction[3], const double closest_pt[3]) const {
    if (gap < 0.0) {
      double contact_force[3] { };
      const double scale = penalty * gap / numNodeFaces;
      for (int i = 0; i < 3; ++i) {
        contact_force[i] = scale * direction[i];
      }
      face.ComputeNodalContactForces(contact_force, closest_pt);
      for (int i = 0; i < 3; ++i) {
        contact_force[i] *= -1.0;
      }
      node.ComputeNodalContactForces(contact_force, closest_pt);
      node.ScatterForceToContactManagerForceVector(contact_manager_force);
      face.ScatterForceToContactManagerForceVector(contact_manager_force);
    }
  }

  double penalty;
  nimble_kokkos::DeviceScalarNodeView contact_manager_force;
};

#endif

class ContactInterface {
 public:
  ContactInterface() = default;
  virtual ~ContactInterface() = default;

  void SetUpPenaltyEnforcement(const double penalty_param) {
#ifdef NIMBLE_HAVE_KOKKOS
    enforcement.penalty = penalty_param;
#endif
  }

#ifdef NIMBLE_HAVE_KOKKOS
  /// \brief Set the contact force manager to a specific Kokkos::View
  void SetContactForce(nimble_kokkos::DeviceScalarNodeView contact_manager_force)
  { enforcement.contact_manager_force = contact_manager_force; }

  void ComputeContact(nimble_kokkos::DeviceContactEntityArrayView contact_nodes,
                      nimble_kokkos::DeviceContactEntityArrayView contact_faces,
                      nimble_kokkos::DeviceScalarNodeView contact_manager_force) {
    ZeroContactForces(contact_manager_force);
    enforcement.contact_manager_force = contact_manager_force;
    DoSearchAndEnforcement(contact_nodes, contact_faces, enforcement);
  }

  inline void ZeroContactForces(nimble_kokkos::DeviceScalarNodeView contact_manager_force) const {
    Kokkos::deep_copy(contact_manager_force, 0.0);
  }

  virtual void DoSearchAndEnforcement(nimble_kokkos::DeviceContactEntityArrayView contact_nodes,
                                      nimble_kokkos::DeviceContactEntityArrayView contact_faces,
                                      PenaltyContactEnforcement enforcement) {
    std::cerr << "Warning: running no-op contact---no interface enabled!" << std::endl;
  }

  KOKKOS_FORCEINLINE_FUNCTION
  void EnforceNodeFaceInteraction(ContactEntity &node, ContactEntity &face, int numNodeFaces, const double gap,
                                  const double direction[3], const double closest_pt[3]) const {
    if (gap < 0.0) {
      enforcement.EnforceContact(node, face, numNodeFaces, gap, direction, closest_pt);
    }
  }

protected:

  PenaltyContactEnforcement enforcement;
#endif
};

}

#endif /* SRC_NIMBLE_CONTACT_INTERFACE_H_ */
