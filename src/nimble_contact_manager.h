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

#include <cfloat>
#include <cmath>
#include <map>
#include <memory>
#include <unordered_map>
#include <stdexcept>
#include <vector>

#include "nimble_contact_entity.h"
#include "nimble_contact_interface.h"
#include "nimble_defs.h"
#include "nimble_exodus_output.h"
#include "nimble_genesis_mesh.h"
#include "nimble_view.h"

#ifdef NIMBLE_HAVE_KOKKOS
#include "nimble_kokkos_contact_defs.h"
#endif

#include "nimble_timer.h"

namespace nimble_kokkos {

class ModelData;

}

namespace nimble {

class DataManager;
class VectorCommunicator;

namespace details {

inline void
getContactForce(const double penalty, const double gap, const double normal[3], double contact_force[3])
{
  const double scale = penalty * gap;
  for (int i = 0; i < 3; ++i) contact_force[i] = scale * normal[i];
}

}  // namespace details

struct PenaltyContactEnforcement
{
  PenaltyContactEnforcement() : penalty(0.0) {}

  template <typename VecType>
  NIMBLE_INLINE_FUNCTION void
  EnforceContact(
      const ContactEntity& node,
      ContactEntity&       face,
      const double         gap,
      const double         normal[3],
      const double         barycentric_coordinates[3],
      VecType&             full_contact_force) const
  {
    constexpr int dim_vec = 3;
    double        contact_force[dim_vec]{};
    details::getContactForce(penalty, gap, normal, contact_force);
    //
    //--- Create a local dummy copy of 'face' to avoid "critical" section
    //
    ContactEntity tmp_face;
    tmp_face.entity_type_                   = face.entity_type_;
    tmp_face.node_id_for_node_1_            = face.node_id_for_node_1_;
    tmp_face.node_id_for_node_2_            = face.node_id_for_node_2_;
    tmp_face.node_id_1_for_fictitious_node_ = face.node_id_1_for_fictitious_node_;
    tmp_face.node_id_2_for_fictitious_node_ = face.node_id_2_for_fictitious_node_;
    tmp_face.node_id_3_for_fictitious_node_ = face.node_id_3_for_fictitious_node_;
    tmp_face.node_id_4_for_fictitious_node_ = face.node_id_4_for_fictitious_node_;
    tmp_face.SetNodalContactForces(contact_force, barycentric_coordinates);
    tmp_face.ScatterForceToContactManagerForceVector<VecType>(full_contact_force);
    //
    //--- Create a local dummy copy of 'node' to avoid "critical" section
    //
    ContactEntity tmp_node;
    tmp_node.entity_type_        = ContactEntity::NODE;
    tmp_node.node_id_for_node_1_ = node.node_id_for_node_1_;
    tmp_node.SetNodalContactForces(contact_force);
    tmp_node.ScatterForceToContactManagerForceVector<VecType>(full_contact_force);
  }

  double penalty;
};

class ContactManager
{
 public:
  typedef enum ProjectionType
  {
    UNKNOWN      = 0,
    NODE_OR_EDGE = 1,
    FACE         = 2
  } PROJECTION_TYPE;

  /// \brief Constructor
  ContactManager(std::shared_ptr<ContactInterface> interface, nimble::DataManager& data_manager);

  /// \brief Default destructor
  virtual ~ContactManager() = default;

  //
  // Interface functions
  //

  /// \brief Indicates whether contact is enabled or not
  ///
  /// \return Boolean for contact
  bool
  ContactEnabled() const
  {
    return contact_enabled_;
  }

  /// \brief Set the penalty parameter for enforcing contact
  ///
  /// \param penalty_parameter Penalty value
  void
  SetPenaltyParameter(double penalty_parameter)
  {
    penalty_parameter_  = penalty_parameter;
    enforcement.penalty = penalty_parameter_;
    contact_interface->SetUpPenaltyEnforcement(penalty_parameter_);
  }

  /// \brief Create contact entities between different blocks
  ///
  /// \param mesh
  /// \param mpi_container
  /// \param primary_block_ids
  /// \param secondary_block_ids
  void
  CreateContactEntities(
      GenesisMesh const&          mesh,
      nimble::VectorCommunicator& mpi_container,
      std::vector<int> const&     primary_block_ids,
      std::vector<int> const&     secondary_block_ids);

  /// \brief Initialize contact visualization output
  ///
  /// \param contact_visualization_exodus_file_name Filename for output
  virtual void
  InitializeContactVisualization(std::string const& contact_visualization_exodus_file_name);

  /// \brief Write information about contact
  ///
  /// \param time_current Current time
  virtual void
  ContactVisualizationWriteStep(double time_current);

  /// \brief Compute the contact force
  ///
  /// \param data_manager
  /// \param step
  /// \param debug_output
  /// \param contact_force
  virtual void
  ComputeContactForce(int step, bool debug_output, nimble::Viewify<2> contact_force);

  /// \brief Returns the number of contact faces
  ///
  /// \return Number of contact faces
  ///
  /// \note When using Kokkos, the data is extracted from the "host".
  size_t
  numContactFaces() const;

  /// \brief Returns the number of contact nodes
  ///
  /// \return Number of contact nodes
  ///
  /// \note When using Kokkos, the data is extracted from the "host".
  size_t
  numContactNodes() const;

  /// \brief Returns the number of contact faces "actively" in collision
  ///
  /// \return Number of active contact faces
  ///
  /// \note When using Kokkos, the data is extracted from the "host".
  size_t
  numActiveContactFaces() const;

  /// \brief Returns the number of contact nodes "actively" in collision
  ///
  /// \return Number of active contact nodes
  ///
  /// \note When using Kokkos, the data is extracted from the "host".
  size_t
  numActiveContactNodes() const;

  /// \brief Return timing information
  /// \return Reference to map of strings to time value
  const std::unordered_map<std::string, double>&
  getTimers();

  //
  // Static functions
  //

  static void
  SkinBlocks(
      GenesisMesh const&             mesh,
      std::vector<int> const&        block_ids,
      int                            entity_id_offset,
      std::vector<std::vector<int>>& skin_faces,
      std::vector<int>&              entity_ids);

  static void
  RemoveInternalSkinFaces(GenesisMesh const& mesh, std::vector<std::vector<int>>& faces, std::vector<int>& entity_ids);

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
  static void
  SimpleClosestPointProjectionSingle(
      const ContactEntity&   node,
      const ContactEntity&   tri,
      PROJECTION_TYPE*       projection_type,
      ContactEntity::vertex* closest_point,
      double&                gap,
      double*                normal,
      double                 tolerance = 1.e-8);

  /// \brief Compute the projection of a point onto a triangular face
  ///
  /// \param[in] node Node to project
  /// \param[in] tri Face to project onto
  /// \param[out] in True if node is inside the facet
  /// \param[out] gap Normal distance when the point projects onto the face
  /// \param[out] normal Unit normal vector outside of face
  /// \param[out] barycentric_coordinates Projection of node on facet
  /// \param[in] tolerance Tolerance to fit into the face (defaut value = 1e-08)
  static void
  Projection(
      const ContactEntity& node,
      const ContactEntity& tri,
      bool&                in,
      double&              gap,
      double*              normal,
      double*              barycentric_coordinates,
      double               tolerance = 1.e-8);

  /// \brief Return the penalty coefficient for enforcing contact force
  ///
  /// \return Penalty value
  double
  GetPenaltyForceParam() const noexcept
  {
    return enforcement.penalty;
  }

 protected:
  //
  // Specific functions
  //

  /// \brief Write visualization data to Exodus file at time t
  ///
  /// \param t Current time
  ///
  /// \note When using Kokkos, the data is extracted from the "host".
  void
  WriteVisualizationData(double t);

  template <typename ArgT>
  void
  CreateContactNodesAndFaces(
      std::vector<std::vector<int>> const& primary_skin_faces,
      std::vector<int> const&              primary_skin_entity_ids,
      std::vector<int> const&              secondary_node_ids,
      std::vector<int> const&              secondary_node_entity_ids,
      std::map<int, double> const&         secondary_node_char_lens,
      ArgT&                                contact_nodes,
      ArgT&                                contact_faces) const;

  void
  BoundingBox(double& x_min, double& x_max, double& y_min, double& y_max, double& z_min, double& z_max) const;

  double
  BoundingBoxAverageCharacteristicLengthOverAllRanks() const;

  void
  ComputeContactForce(int step, bool debug_output)
  {
    if (penalty_parameter_ <= 0.0) {
      throw std::invalid_argument("\nError in ComputeContactForce(), invalid penalty_parameter.\n");
    }
    double background_grid_cell_size = BoundingBoxAverageCharacteristicLengthOverAllRanks();
#ifdef NIMBLE_HAVE_KOKKOS
    contact_interface->ComputeContact(contact_nodes_d_, contact_faces_d_, force_d_);
#endif
  }

  void
  ApplyDisplacements(nimble_kokkos::DeviceVectorNodeView displacement_d);
  void
  GetForces(nimble_kokkos::DeviceVectorNodeView contact_force_d) const;

  void
  BruteForceBoxIntersectionSearch(std::vector<ContactEntity> const& nodes, std::vector<ContactEntity> const& triangles);

  void
  ClosestPointProjection(
      const ContactEntity*   nodes,
      const ContactEntity*   triangles,
      ContactEntity::vertex* closest_points,
      PROJECTION_TYPE*       projection_types,
      std::size_t            num_elements);

  void
  ClosestPointProjectionSingle(
      const ContactEntity&   node,
      const ContactEntity&   tri,
      ContactEntity::vertex* closest_point,
      PROJECTION_TYPE*       projection_type,
      double                 tolerance);

  /// \brief Returns a read-only reference to contact face entity
  ///
  /// \param i_face Index of the contact face
  /// \return Read-only reference to contact face entity
  ///
  /// \note When using Kokkos, the data is extracted from the "host".
  const ContactEntity&
  getContactFace(size_t i_face) const;

  /// \brief Returns a read-only reference to contact node entity
  ///
  /// \param i_node Index of the contact node
  /// \return Read-only reference to contact node entity
  ///
  /// \note When using Kokkos, the data is extracted from the "host".
  const ContactEntity&
  getContactNode(size_t i_node) const;

  /// \brief Reset the contact status for all entities
  void
  ResetContactStatus();

  /// \brief Zero the contact forces
  void
  ZeroContactForce();

  template <typename VecType>
  void
  EnforceNodeFaceInteraction(
      ContactEntity& node,
      ContactEntity& face,
      const double   gap,
      const double   direction[3],
      const double   facet_coordinates[3],
      VecType&       full_contact_force) const
  {
    enforcement.EnforceContact<VecType>(node, face, gap, direction, facet_coordinates, full_contact_force);
  }

  NIMBLE_INLINE_FUNCTION
  void
  startTimer(const std::string& str_val)
  {
#ifdef NIMBLE_TIME_CONTACT
    watch_.Start(str_val);
#endif
  }

  NIMBLE_INLINE_FUNCTION
  void
  stopTimer(const std::string& str_val)
  {
#ifdef NIMBLE_TIME_CONTACT
    watch_.Stop(str_val);
#endif
  }
  //--- Variables

  nimble::DataManager& data_manager_;

  nimble::TimeKeeper total_search_time;
  nimble::TimeKeeper total_enforcement_time;
  std::size_t        total_num_contacts = 0;

  PenaltyContactEnforcement enforcement;

  bool   contact_enabled_ = false;
  double penalty_parameter_;

  std::vector<int>           node_ids_;
  std::vector<double>        model_coord_;
  std::vector<double>        coord_;
  std::vector<double>        force_;
  std::vector<ContactEntity> contact_faces_;
  std::vector<ContactEntity> contact_nodes_;

  double               contact_visualization_model_coord_bounding_box_[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  nimble::GenesisMesh  genesis_mesh_for_contact_visualization_;
  nimble::ExodusOutput exodus_output_for_contact_visualization_;

  nimble_kokkos::DeviceIntegerArrayView node_ids_d_  = nimble_kokkos::DeviceIntegerArrayView("contact node_ids_d", 1);
  nimble_kokkos::DeviceScalarNodeView model_coord_d_ = nimble_kokkos::DeviceScalarNodeView("contact model_coord_d", 1);
  nimble_kokkos::DeviceScalarNodeView coord_d_       = nimble_kokkos::DeviceScalarNodeView("contact coord_d", 1);
  nimble_kokkos::DeviceScalarNodeView force_d_       = nimble_kokkos::DeviceScalarNodeView("contact force_d", 1);

  nimble_kokkos::DeviceContactEntityArrayView contact_faces_d_ =
      nimble_kokkos::DeviceContactEntityArrayView("contact_faces_d", 1);
  nimble_kokkos::DeviceContactEntityArrayView contact_nodes_d_ =
      nimble_kokkos::DeviceContactEntityArrayView("contact_nodes_d", 1);

  // TODO remove this once enforcement is on device
  nimble_kokkos::HostContactEntityArrayView contact_faces_h_ =
      nimble_kokkos::HostContactEntityArrayView("contact_faces_h", 1);
  nimble_kokkos::HostContactEntityArrayView contact_nodes_h_ =
      nimble_kokkos::HostContactEntityArrayView("contact_nodes_h", 1);

  nimble_kokkos::ModelData* model_data_ = nullptr;

  std::unordered_map<std::string, double> timers_;
#ifdef NIMBLE_TIME_CONTACT
  nimble::Timer watch_;
#endif

  std::shared_ptr<ContactInterface> contact_interface;
};

std::shared_ptr<nimble::ContactManager>
GetContactManager(std::shared_ptr<ContactInterface> interface, nimble::DataManager& data_manager);

}  // namespace nimble

#endif  // NIMBLE_MATERIAL_H
