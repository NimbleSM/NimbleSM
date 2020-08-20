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
#include <cfloat>
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

namespace details {

inline
void getContactForce(
    const double penalty,
    const double gap,
    const double normal[3],
    double contact_force[3]
)
{
  const double scale = penalty * gap;
  for (int i = 0; i < 3; ++i)
    contact_force[i] = scale * normal[i];
}

}

struct PenaltyContactEnforcement {
  PenaltyContactEnforcement() : penalty(0.0) {}

  template < typename VecType >
  NIMBLE_INLINE_FUNCTION
  void EnforceContact(ContactEntity &node, 
                      ContactEntity &face, 
                      const double gap,
                      const double normal[3], 
                      const double barycentric_coordinates[3],
                      VecType &full_contact_force) const {
      double contact_force[3] { };
      details::getContactForce(penalty, gap, normal, contact_force);
      //
      face.SetNodalContactForces(contact_force, barycentric_coordinates);
      node.SetNodalContactForces(contact_force);
      //
      node.ScatterForceToContactManagerForceVector<VecType>(full_contact_force);
      face.ScatterForceToContactManagerForceVector<VecType>(full_contact_force);
  }

  double penalty;

};

class ContactManager {

public:

    typedef enum ProjectionType {
      UNKNOWN=0,
      NODE_OR_EDGE=1,
      FACE=2
    } PROJECTION_TYPE;

    explicit ContactManager( std::shared_ptr<ContactInterface> interface );

    virtual ~ContactManager() = default;

    bool ContactEnabled() const { return contact_enabled_; }

    static void SkinBlocks(GenesisMesh const & mesh,
                    std::vector<int> const & block_ids,
                    int entity_id_offset,
                    std::vector< std::vector<int> > & skin_faces,
                    std::vector<int> & entity_ids);

    static void RemoveInternalSkinFaces(GenesisMesh const & mesh,
                                 std::vector< std::vector<int> >& faces,
                                 std::vector<int>& entity_ids);

    void SetPenaltyParameter(double penalty_parameter) {
      penalty_parameter_ = penalty_parameter;
      enforcement.penalty = penalty_parameter_;
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

    void ApplyDisplacements(const double * displacement);
    void GetForces(double * contact_force) const;

#ifdef NIMBLE_HAVE_KOKKOS
    // Kokkos versions of ApplyDisplacements and GetForces
    void ApplyDisplacements(nimble_kokkos::DeviceVectorNodeView displacement_d);
    void GetForces(nimble_kokkos::DeviceVectorNodeView contact_force_d) const;
#endif

    virtual void ComputeContactForce(int step, bool debug_output);

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


// DEPRECATED
    /// \brief Compute the projection of a point onto a triangular face
    ///
    /// \param[in] node Node to project
    /// \param[in] tri Face to project onto
    /// \param[out] projection_type Result of projection (FACE: success,
    ///             UNKNOWN: point projects outside the face)
    /// \param[out] closest_point Projection
    /// \param[out] gap Normal distance when the point projects onto the face
    /// \param[out] normal Unit normal vector outside of face
    /// \param[in] tolerance Tolerance to fit into the face (defaut value = 1e-08)
    static void SimpleClosestPointProjectionSingle( const ContactEntity &node,
        const ContactEntity &tri,
        PROJECTION_TYPE *projection_type,
        ContactEntity::vertex *closest_point,
        double &gap,
        double *normal,
        double tolerance = 1.e-8);

    /// \brief Compute the projection of a point onto a triangular face
    ///
    /// \param[in] node Node to project
    /// \param[in] tri Face to project onto
    /// \param[out] in True if node is inside the facet 
    /// \param[out] gap Normal distance when the point projects onto the face
    /// \param[out] normal Unit normal vector outside of face
    /// \param[out] barycentric_coordinates Projection of node on facet
    /// \param[in] tolerance Tolerance to fit into the face (defaut value = 1e-08)
    static void Projection( const ContactEntity &node,
        const ContactEntity &tri,
        bool   &in,
        double &gap,
        double *normal,
        double *barycentric_coordinates,
        double tolerance = 1.e-8);

    virtual void InitializeContactVisualization(std::string const & contact_visualization_exodus_file_name);

    virtual void ContactVisualizationWriteStep(double time_current);

    /// \brief Returns the number of contact faces "actively" in collision
    ///
    /// \return Number of active contact faces
    ///
    /// \note When using Kokkos, the data is extracted from the "host".
    std::size_t numActiveContactFaces() const {
       std::size_t num_contacts = 0;
       for (size_t i = 0; i < numContactFaces(); ++i) {
         const auto &myface = getContactFace(i);
         num_contacts += static_cast<size_t>(myface.contact_status());
       }
       return num_contacts;
    }

    /// \brief Returns a read-only reference to contact face entity
    ///
    /// \param i_face Index of the contact face
    /// \return Read-only reference to contact face entity
    ///
    /// \note When using Kokkos, the data is extracted from the "host".
    const ContactEntity& getContactFace(size_t i_face) const {
#ifdef NIMBLE_HAVE_KOKKOS
      return contact_faces_h_(i_face);
#else
      return contact_faces_[i_face];
#endif
    }

    /// \brief Returns the number of contact faces
    ///
    /// \return Number of contact faces
    ///
    /// \note When using Kokkos, the data is extracted from the "host".
    size_t numContactFaces() const {
#ifdef NIMBLE_HAVE_KOKKOS
    return contact_faces_h_.extent(0);
#else
    return contact_faces_.size();
#endif
    }

    /// \brief Returns a read-only reference to contact node entity
    ///
    /// \param i_node Index of the contact node
    /// \return Read-only reference to contact node entity
    ///
    /// \note When using Kokkos, the data is extracted from the "host".
    const ContactEntity& getContactNode(size_t i_node) const {
#ifdef NIMBLE_HAVE_KOKKOS
      return contact_nodes_h_(i_node);
#else
      return contact_nodes_[i_node];
#endif
    }

    /// \brief Returns the number of contact nodes
    ///
    /// \return Number of contact nodes
    ///
    /// \note When using Kokkos, the data is extracted from the "host".
    size_t numContactNodes() const {
#ifdef NIMBLE_HAVE_KOKKOS
      return contact_nodes_h_.extent(0);
#else
      return contact_nodes_.size();
#endif
    }

    /// \brief Zero the contact forces
    void zeroContactForce();

protected:
 
  PenaltyContactEnforcement enforcement;

  template< typename VecType>
  void EnforceNodeFaceInteraction(
      ContactEntity &node, 
      ContactEntity &face, 
      const double gap,
      const double direction[3], 
      const double facet_coordinates[3],
      VecType &full_contact_force) const {
    enforcement.EnforceContact<VecType>(node, face, gap, direction, facet_coordinates, 
                full_contact_force);
  }

  /// Routine to write the contact data to Exodus file at time t
  ///
  /// \param t Time for the current data
  ///
  /// \note When using Kokkos, the data is extracted from the "host".
  void WriteVisualizationData(double t);

  //--- Variables

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
