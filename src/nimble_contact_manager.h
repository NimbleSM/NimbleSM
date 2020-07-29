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

#ifndef NIMBLE_CONTACT_H
#define NIMBLE_CONTACT_H

#include <vector>
#include <map>
#include <memory>
#include <float.h>
#include <cmath>

#include "nimble_contact_entity.h"
#include "nimble_contact_interface.h"
#include "nimble_genesis_mesh.h"
#include "nimble_exodus_output.h"
#include "nimble.mpi.utils.h"

#ifdef NIMBLE_HAVE_KOKKOS
  #include "nimble_kokkos_defs.h"
  #include "nimble_kokkos_contact_defs.h"
#endif

#ifdef NIMBLE_HAVE_BVH
  #include <bvh/kdop.hpp>
  #include <bvh/tree.hpp>
  #include <bvh/patch.hpp>
  #include <bvh/perf/instrument.hpp>
#ifdef BVH_ENABLE_VT
  #include <bvh/vt/collection.hpp>
  #include <bvh/vt/helpers.hpp>
  #include <bvh/vt/collision_world.hpp>
#endif
#endif

namespace nimble {

class ContactInterface;

  class ContactManager {

  public:

    typedef enum ProjectionType {
      UNKNOWN=0,
      NODE_OR_EDGE=1,
      FACE=2
    } PROJECTION_TYPE;

    explicit ContactManager( std::shared_ptr<ContactInterface> interface );

    virtual ~ContactManager() {}

    bool ContactEnabled() { return contact_enabled_; }

    void SkinBlocks(GenesisMesh const & mesh,
                    std::vector<int> const & block_ids,
                    int entity_id_offset,
                    std::vector< std::vector<int> > & skin_faces,
                    std::vector<int> & entity_ids);

    void RemoveInternalSkinFaces(GenesisMesh const & mesh,
                                 std::vector< std::vector<int> >& faces,
                                 std::vector<int>& entity_ids);

    void SetPenaltyParameter(double penalty_parameter) {
      penalty_parameter_ = penalty_parameter;
      contact_interface->SetUpPenaltyEnforcement(penalty_parameter_);
    }

    void CreateContactEntities(GenesisMesh const & mesh,
                               std::vector<int> const & master_block_ids,
                               std::vector<int> const & slave_block_ids) {
      nimble::MPIContainer mpi_container;
      CreateContactEntities(mesh, mpi_container, master_block_ids, slave_block_ids);
    }


    void CreateContactEntities(GenesisMesh const & mesh,
                               nimble::MPIContainer & mpi_container,
                               std::vector<int> const & master_block_ids,
                               std::vector<int> const & slave_block_ids);

    template <typename ArgT>
    void CreateContactNodesAndFaces(std::vector< std::vector<int> > const & master_skin_faces,
                                    std::vector<int> const & master_skin_entity_ids,
                                    std::vector<int> const & slave_node_ids,
                                    std::vector<int> const & slave_skin_entity_ids,
                                    std::map<int, double> const & slave_node_char_lens,
                                    ArgT& contact_nodes,
                                    ArgT& contact_faces) const ;

    void BoundingBox(double& x_min,
                     double& x_max,
                     double& y_min,
                     double& y_max,
                     double& z_min,
                     double& z_max) const ;

    double BoundingBoxAverageCharacteristicLengthOverAllRanks() const ;

    void ApplyDisplacements(const double * const displacement) {
      for (unsigned int i_node=0; i_node<node_ids_.size() ; i_node++) {
        int node_id = node_ids_[i_node];
        for (int i=0 ; i<3 ; i++) {
          coord_[3*i_node + i] = model_coord_[3*i_node + i] + displacement[3*node_id + i];
        }
      }
      for (unsigned int i_face=0 ; i_face<contact_faces_.size() ; i_face++) {
        contact_faces_[i_face].SetCoordinates(coord_.data());
      }
      for (unsigned int i_node=0 ; i_node<contact_nodes_.size() ; i_node++) {
        contact_nodes_[i_node].SetCoordinates(coord_.data());
      }
    }

#ifdef NIMBLE_HAVE_KOKKOS

    void ApplyDisplacements(nimble_kokkos::DeviceVectorNodeView displacement_d) {

      int num_nodes_in_contact_manager = node_ids_d_.extent(0);
      int num_contact_node_entities = contact_nodes_d_.extent(0);
      int num_contact_face_entities = contact_faces_d_.extent(0);

      // circumvent lambda *this glitch
      nimble_kokkos::DeviceIntegerArrayView node_ids = node_ids_d_;
      nimble_kokkos::DeviceScalarNodeView model_coord = model_coord_d_;
      nimble_kokkos::DeviceScalarNodeView coord = coord_d_;
      nimble_kokkos::DeviceContactEntityArrayView contact_nodes = contact_nodes_d_;
      nimble_kokkos::DeviceContactEntityArrayView contact_faces = contact_faces_d_;

      Kokkos::parallel_for("ContactManager::ApplyDisplacements set coord_d_ vector",
                           num_nodes_in_contact_manager,
                           KOKKOS_LAMBDA(const int i) {
        int node_id = node_ids(i);
        coord(3*i)   = model_coord(3*i)   + displacement_d(node_id, 0);
        coord(3*i+1) = model_coord(3*i+1) + displacement_d(node_id, 1);
        coord(3*i+2) = model_coord(3*i+2) + displacement_d(node_id, 2);
      });

      Kokkos::parallel_for("ContactManager::ApplyDisplacements set contact node entity displacements",
                         num_contact_node_entities,
                         KOKKOS_LAMBDA(const int i_node) {
        contact_nodes(i_node).SetCoordinates(coord);
      });

      Kokkos::parallel_for("ContactManager::ApplyDisplacements set contact face entity displacements",
                         num_contact_face_entities,
                         KOKKOS_LAMBDA(const int i_face) {
        contact_faces(i_face).SetCoordinates(coord);
      });

    }
#endif

    void GetForces(double * const contact_force) {
      for (unsigned int i_node=0; i_node<node_ids_.size() ; i_node++) {
        int node_id = node_ids_[i_node];
        for (int i=0 ; i<3 ; i++) {
          contact_force[3*node_id + i] = force_[3*i_node + i];
        }
      }
    }

#ifdef NIMBLE_HAVE_KOKKOS
    void GetForces(nimble_kokkos::DeviceVectorNodeView contact_force_d) {

      int num_nodes_in_contact_manager = node_ids_d_.extent(0);

      // circumvent lambda *this glitch
      nimble_kokkos::DeviceIntegerArrayView node_ids = node_ids_d_;
      nimble_kokkos::DeviceScalarNodeView force = force_d_;

      Kokkos::parallel_for("ContactManager::GetForces",
                         num_nodes_in_contact_manager,
                         KOKKOS_LAMBDA(const int i) {
        int node_id = node_ids(i);
        contact_force_d(node_id, 0) = force(3*i);
        contact_force_d(node_id, 1) = force(3*i+1);
        contact_force_d(node_id, 2) = force(3*i+2);
      });
    }
#endif

    void ComputeContactForce(int step, bool debug_output);

    void BruteForceBoxIntersectionSearch(std::vector<ContactEntity> const & nodes,
                                         std::vector<ContactEntity> const & triangles);

    void ClosestPointProjection(const ContactEntity *nodes,
                                const ContactEntity *triangles,
                                ContactEntity::vertex *closest_points,
                                PROJECTION_TYPE *projection_types,
                                std::size_t num_elements);
    void ClosestPointProjectionSingle( const ContactEntity &node,
        const ContactEntity &tri,
        ContactEntity::vertex *closest_point,
        PROJECTION_TYPE *projection_type,
        double tolerance);
    void SimpleClosestPointProjectionSingle( const ContactEntity &node,
        const ContactEntity &tri,
        ContactEntity::vertex *closest_point,
        PROJECTION_TYPE *projection_type,
        double tolerance);

    virtual void InitializeContactVisualization(std::string const & contact_visualization_exodus_file_name);

    virtual void ContactVisualizationWriteStep(double time_current);

  protected:

    void InitializeContactVisualizationImpl(std::string const & contact_visualization_exodus_file_name,
                                            nimble::GenesisMesh &mesh,
                                            nimble::ExodusOutput &out,
                                            ContactEntity *faces, std::size_t nfaces,
                                            ContactEntity *nodes, std::size_t nnodes);
    void WriteVisualizationData( double t, nimble::GenesisMesh &mesh,
                                 nimble::ExodusOutput &out,
        ContactEntity *faces, std::size_t nfaces,
        ContactEntity *nodes, std::size_t nnodes );

    bool contact_enabled_ = false;
    double penalty_parameter_;

    std::vector<int> node_ids_;
    std::vector<double> model_coord_;
    std::vector<double> coord_;
    std::vector<double> force_;
    std::vector<ContactEntity> contact_faces_;
    std::vector<ContactEntity> contact_nodes_;

    double contact_visualization_model_coord_bounding_box_[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    nimble::GenesisMesh genesis_mesh_for_contact_visualization_;
    nimble::ExodusOutput exodus_output_for_contact_visualization_;

#ifdef NIMBLE_HAVE_KOKKOS
    nimble_kokkos::DeviceIntegerArrayView node_ids_d_ = nimble_kokkos::DeviceIntegerArrayView("contact node_ids_d", 1);
    nimble_kokkos::DeviceScalarNodeView model_coord_d_ = nimble_kokkos::DeviceScalarNodeView("contact model_coord_d", 1);
    nimble_kokkos::DeviceScalarNodeView coord_d_ = nimble_kokkos::DeviceScalarNodeView("contact coord_d", 1);
    nimble_kokkos::DeviceScalarNodeView force_d_ = nimble_kokkos::DeviceScalarNodeView("contact force_d", 1);

    nimble_kokkos::DeviceContactEntityArrayView contact_faces_d_ = nimble_kokkos::DeviceContactEntityArrayView("contact_faces_d", 1);
    nimble_kokkos::DeviceContactEntityArrayView contact_nodes_d_ = nimble_kokkos::DeviceContactEntityArrayView("contact_nodes_d", 1);

    // TODO remove this once enforcement is on device
    nimble_kokkos::HostContactEntityArrayView contact_faces_h_ = nimble_kokkos::HostContactEntityArrayView("contact_faces_h", 1);
    nimble_kokkos::HostContactEntityArrayView contact_nodes_h_ = nimble_kokkos::HostContactEntityArrayView("contact_nodes_h", 1);
#endif

    std::shared_ptr<ContactInterface> contact_interface;
};

} // namespace nimble

#endif // NIMBLE_MATERIAL_H
