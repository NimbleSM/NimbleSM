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
#include <float.h>
#include <cmath>
#include "nimble_genesis_mesh.h"
#include "nimble_kokkos_defs.h"

#ifdef NIMBLE_HAVE_EXTRAS
  #include "nimble_extras_contact_includes.h"
#endif

#ifdef NIMBLE_HAVE_BVH
  #include <bvh/kdop.hpp>
  #include <bvh/tree.hpp>
#endif

namespace nimble {

  NIMBLE_INLINE_FUNCTION
  double TriangleArea(const double * const pt_1,
                      const double * const pt_2,
                      const double * const pt_3) {

    double a[3], b[3], cross[3], area;
    for (int i=0 ; i<3 ; ++i) {
      a[i] = pt_2[i] - pt_1[i];
      b[i] = pt_3[i] - pt_1[i];
    }
    cross[0] = b[1]*a[2] - b[2]*a[1];
    cross[1] = b[0]*a[2] - b[2]*a[0];
    cross[2] = b[0]*a[1] - b[1]*a[0];
    area = 0.5 * std::sqrt(cross[0]*cross[0] + cross[1]*cross[1] + cross[2]*cross[2]);
    return area;
  }

  void ParseContactCommand(std::string const & command,
                           std::vector<std::string> & master_block_names,
                           std::vector<std::string> & slave_block_names,
                           double & penalty_parameter);

  class ContactEntity {

  public:

    typedef enum ContactEntityType {
      NONE=0,
      NODE=1,
      TRIANGLE=3
    } CONTACT_ENTITY_TYPE;

    ContactEntity() : entity_type_(NONE) {}

    ContactEntity(CONTACT_ENTITY_TYPE entity_type,
                  int contant_entity_global_id,
                  double const coord[],
                  double characteristic_length,
                  int node_id_for_node_1,
                  int node_id_for_node_2 = 0,
                  int node_ids_for_fictitious_node[4] = 0)
      : entity_type_(entity_type),
        contant_entity_global_id_(contant_entity_global_id),
        char_len_(characteristic_length) {

          // contact entities must be either nodes (one node) or trianglular faces (three nodes)
          if (entity_type_ == NODE) {
            num_nodes_ = 1;
            node_id_for_node_1_ = node_id_for_node_1;
          }
          else if (entity_type_ == TRIANGLE) {
            num_nodes_ = 3;
            node_id_for_node_1_ = node_id_for_node_1;
            node_id_for_node_2_ = node_id_for_node_2;
            for (int i=0 ; i<4 ; i++) {
              node_ids_for_fictitious_node_[i] = node_ids_for_fictitious_node[i];
            }
          }
          else{
            printf("\n**** Error in ContactEntity constructor, invalid entity type.\n");
          }

          for (unsigned int i=0 ; i<3*num_nodes_ ; i++) {
            coord_[i] = coord[i];
          }

          SetBoundingBox();
        }

    virtual ~ContactEntity() {}

    template <typename ArgT>
    void SetCoordinates(ArgT coord) {
      int n0 = 3*node_id_for_node_1_;
      coord_[0] = coord[n0];
      coord_[1] = coord[n0+1];
      coord_[2] = coord[n0+2];
      if (entity_type_ == TRIANGLE) {
        n0 = 3*node_id_for_node_2_;
        coord_[3] = coord[n0];
        coord_[4] = coord[n0+1];
        coord_[5] = coord[n0+2];
        n0 = 3*node_ids_for_fictitious_node_[0];
        int n1 = 3*node_ids_for_fictitious_node_[1];
        int n2 = 3*node_ids_for_fictitious_node_[2];
        int n3 = 3*node_ids_for_fictitious_node_[3];
        coord_[6] = (coord[n0]   + coord[n1]   + coord[n2]   + coord[n3]  ) / 4.0;
        coord_[7] = (coord[n0+1] + coord[n1+1] + coord[n2+1] + coord[n3+1]) / 4.0;
        coord_[8] = (coord[n0+2] + coord[n1+2] + coord[n2+2] + coord[n3+2]) / 4.0;
      }
      SetBoundingBox();
    }

    template <typename ArgT>
    void GetForces(ArgT force) {
      int n = 3*node_id_for_node_1_;
      force[n]   += force_[0];
      force[n+1] += force_[1];
      force[n+2] += force_[2];
      if (entity_type_ == TRIANGLE) {
        n = 3*node_id_for_node_2_;
        force[n]   += force_[3];
        force[n+1] += force_[4];
        force[n+2] += force_[5];
        for (int i=0 ; i<4 ; i++) {
          int n = 3*node_ids_for_fictitious_node_[i];
          force[n]   += force_[6] / 4.0;
          force[n+1] += force_[7] / 4.0;
          force[n+2] += force_[8] / 4.0;
        }
      }
    }

    // functions required for NimbleSMExtras contact search
    double get_x_min() const { return bounding_box_x_min_; }
    double get_x_max() const { return bounding_box_x_max_; }
    double get_y_min() const { return bounding_box_y_min_; }
    double get_y_max() const { return bounding_box_y_max_; }
    double get_z_min() const { return bounding_box_z_min_; }
    double get_z_max() const { return bounding_box_z_max_; }

    // Functions for bvh contact search
    int contant_entity_global_id() const { return contant_entity_global_id_; }
    const double* centroid() const { return centroid_; }
#ifdef NIMBLE_HAVE_BVH
    bvh::dop_26<double> kdop() const {
      std::vector< std::array< double, 3 > > vertices;
      for ( std::size_t i = 0; i < 3*num_nodes_; i += 3 )
      {
        vertices.push_back( {{ coord_[i], coord_[i + 1], coord_[i + 2] }} );
      }
      const double inflation_length = 0.15 * char_len_;
      return bvh::dop_26< double >::from_vertices( vertices.begin(), vertices.end(), inflation_length );
    }
#endif

    void SetBoundingBox() {

      centroid_[0] = centroid_[1] = centroid_[2] = 0.0;

      bounding_box_x_min_ = DBL_MAX;
      bounding_box_x_max_ = -DBL_MAX;
      bounding_box_y_min_ = DBL_MAX;
      bounding_box_y_max_ = -DBL_MAX;
      bounding_box_z_min_ = DBL_MAX;
      bounding_box_z_max_ = -DBL_MAX;

      for (int i_node=0 ; i_node<num_nodes_ ; i_node++) {

        centroid_[0] += coord_[3*i_node];
        centroid_[1] += coord_[3*i_node+1];
        centroid_[2] += coord_[3*i_node+2];

        if (coord_[3*i_node] < bounding_box_x_min_)
          bounding_box_x_min_ = coord_[3*i_node];
        if (coord_[3*i_node] > bounding_box_x_max_)
          bounding_box_x_max_ = coord_[3*i_node];
        if (coord_[3*i_node+1] < bounding_box_y_min_)
          bounding_box_y_min_ = coord_[3*i_node+1];
        if (coord_[3*i_node+1] > bounding_box_y_max_)
          bounding_box_y_max_ = coord_[3*i_node+1];
        if (coord_[3*i_node+2] < bounding_box_z_min_)
          bounding_box_z_min_ = coord_[3*i_node+2];
        if (coord_[3*i_node+2] > bounding_box_z_max_)
          bounding_box_z_max_ = coord_[3*i_node+2];
      }

      centroid_[0] /= num_nodes_;
      centroid_[1] /= num_nodes_;
      centroid_[2] /= num_nodes_;

      double inflation_length = 0.15 * char_len_;

      bounding_box_x_min_ -= inflation_length;
      bounding_box_x_max_ += inflation_length;
      bounding_box_y_min_ -= inflation_length;
      bounding_box_y_max_ += inflation_length;
      bounding_box_z_min_ -= inflation_length;
      bounding_box_z_max_ += inflation_length;
    }

    void ComputeNodalContactForces(const double * const contact_force,
                                   const double * const closest_point_projection) {

      double N[3] = {0.0, 0.0, 0.0};

      if (entity_type_ == NODE) {
        N[0] = 1.0;
      }
      else if (entity_type_ == TRIANGLE) {
        // Find the natural coordinates of the closest_point_projection
        double area_1 = TriangleArea(closest_point_projection, &coord_[3], &coord_[6]);
        double area_2 = TriangleArea(&coord_[0], closest_point_projection, &coord_[6]);
        double full_area = TriangleArea(&coord_[0], &coord_[3], &coord_[6]);
        N[0] = area_1/full_area;
        N[1] = area_2/full_area;
        N[2] = 1.0 - N[0] - N[1];
      }

      for (int i_node = 0; i_node < num_nodes_ ; ++i_node) {
        for (int i=0 ; i<3 ; ++i) {
          force_[3*i_node + i] = N[i_node] * contact_force[i];
        }
      }
    }

    int contant_entity_global_id_ = -1;

    // positions of nodes that define triangular contact patch
    CONTACT_ENTITY_TYPE entity_type_ = NONE;
    int num_nodes_ = 0;
    double coord_[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double force_[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double char_len_ = 0.0;

    // bounding box for NimbleSMExtras contact search routines
    double bounding_box_x_min_ = -DBL_MAX;
    double bounding_box_x_max_ =  DBL_MAX;
    double bounding_box_y_min_ = -DBL_MAX;
    double bounding_box_y_max_ =  DBL_MAX;
    double bounding_box_z_min_ = -DBL_MAX;
    double bounding_box_z_max_ =  DBL_MAX;

    // centroid for bvh search routines
    double centroid_[3] = {0.0, 0.0, 0.0};

    // map for moving displacement/forces to/from the contact manager data structures
    int node_id_for_node_1_ = -1;
    int node_id_for_node_2_ = -1;
    int node_ids_for_fictitious_node_[4] = {-1, -1, -1, -1}; // fictitious node maps to multiple real nodes
  };

  class ContactManager {

  public:

    ContactManager() : penalty_parameter_(0.0)
#ifdef NIMBLE_HAVE_EXTRAS
      , contact_nodes_search_tree_("contact nodes search tree")
      , contact_faces_search_tree_("contact faces search tree")
#endif
      {}

    virtual ~ContactManager() {}

    bool ContactEnabled() { return contact_enabled_; }

    std::vector< std::vector<int> > SkinBlocks(GenesisMesh const & mesh,
                                               std::vector<int> block_ids);

    void SetPenaltyParameter(double penalty_parameter) {
      penalty_parameter_ = penalty_parameter;
    }

    void CreateContactEntities(GenesisMesh const & mesh,
                               std::vector<int> master_block_ids,
                               std::vector<int> slave_block_ids);

    template <typename ArgT>
    void CreateContactNodesAndFaces(std::vector< std::vector<int> > const & master_skin_faces,
                                    std::vector<int> const & slave_node_ids,
                                    std::map<int, double> const & slave_node_char_lens,
                                    ArgT& contact_nodes,
                                    ArgT& contact_faces) const ;

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

      Kokkos::parallel_for("ContactManager::ApplyDisplacements set coord_d_ vector",
                         Kokkos::RangePolicy<nimble_kokkos::kokkos_device_execution_space>(0, node_ids_d_.extent(0)),
                         KOKKOS_LAMBDA(const int i) {
        int node_id = node_ids_d_(i);
        coord_d_(3*i)   = model_coord_d_(3*i)   + displacement_d(node_id, 0);
        coord_d_(3*i+1) = model_coord_d_(3*i+1) + displacement_d(node_id, 1);
        coord_d_(3*i+2) = model_coord_d_(3*i+2) + displacement_d(node_id, 2);
      });
      Kokkos::fence();
      Kokkos::parallel_for("ContactManager::ApplyDisplacements set contact node entity displacements",
                         Kokkos::RangePolicy<nimble_kokkos::kokkos_device_execution_space>(0, contact_nodes_d_.extent(0)),
                         KOKKOS_LAMBDA(const int i_node) {
        contact_nodes_d_(i_node).SetCoordinates(coord_d_);
      });

      Kokkos::parallel_for("ContactManager::ApplyDisplacements set contact face entity displacements",
                         Kokkos::RangePolicy<nimble_kokkos::kokkos_device_execution_space>(0, contact_faces_d_.extent(0)),
                         KOKKOS_LAMBDA(const int i_face) {
        contact_faces_d_(i_face).SetCoordinates(coord_d_);
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

      Kokkos::parallel_for("ContactManager::GetForces",
                         Kokkos::RangePolicy<nimble_kokkos::kokkos_device_execution_space>(0, node_ids_d_.extent(0)),
                         KOKKOS_LAMBDA(const int i) {
        int node_id = node_ids_d_(i);
        contact_force_d(node_id, 0) = force_d_(3*i);
        contact_force_d(node_id, 1) = force_d_(3*i+1);
        contact_force_d(node_id, 2) = force_d_(3*i+2);
      });
    }
#endif

    void ComputeContactForce(int step, bool debug_output);

    void WriteContactEntitiesToVTKFile(int step);
    void WriteContactEntitiesToVTKFile(const std::vector<ContactEntity> &_faces,
                                       const std::vector<ContactEntity> &_nodes,
                                       const std::string &_prefix, int step);

#ifdef NIMBLE_HAVE_BVH
    void VisualizeCollisionInfo(const bvh::bvh_tree_26d &faces_tree, const bvh::bvh_tree_26d &nodes_tree,
                                const bvh::bvh_tree_26d::collision_query_result_type &collision_result,
                                int step);
#endif

  protected:

    bool contact_enabled_ = false;
    double penalty_parameter_;

    std::vector<int> node_ids_;
    std::vector<double> model_coord_;
    std::vector<double> coord_;
    std::vector<double> force_;
    std::vector<ContactEntity> contact_faces_;
    std::vector<ContactEntity> contact_nodes_;

#ifdef NIMBLE_HAVE_KOKKOS
    nimble_kokkos::DeviceIntegerArrayView node_ids_d_ = nimble_kokkos::DeviceIntegerArrayView("contact node_ids_d", 1);
    nimble_kokkos::DeviceScalarNodeView model_coord_d_ = nimble_kokkos::DeviceScalarNodeView("contact model_coord_d", 1);
    nimble_kokkos::DeviceScalarNodeView coord_d_ = nimble_kokkos::DeviceScalarNodeView("contact coord_d", 1);
    nimble_kokkos::DeviceScalarNodeView force_d_ = nimble_kokkos::DeviceScalarNodeView("contact force_d", 1);

    using HostContactEntityArrayView = Kokkos::View< ContactEntity*, nimble_kokkos::kokkos_layout, nimble_kokkos::kokkos_host >;
    using DeviceContactEntityArrayView = Kokkos::View< ContactEntity*, nimble_kokkos::kokkos_layout, nimble_kokkos::kokkos_device >;
    DeviceContactEntityArrayView contact_faces_d_ = DeviceContactEntityArrayView("contact_faces_d", 1);
    DeviceContactEntityArrayView contact_nodes_d_ = DeviceContactEntityArrayView("contact_nodes_d", 1);
#endif

#ifdef NIMBLE_HAVE_EXTRAS
    stk::search::mas_aabb_tree<double, nimble_kokkos::kokkos_device_execution_space> contact_nodes_search_tree_;
    stk::search::mas_aabb_tree<double, nimble_kokkos::kokkos_device_execution_space> contact_faces_search_tree_;
#endif
  };

} // namespace nimble

#endif // NIMBLE_MATERIAL_H
