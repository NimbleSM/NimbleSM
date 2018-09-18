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
  double TriangleArea(double pt_1_x, double pt_1_y, double pt_1_z,
                      double pt_2_x, double pt_2_y, double pt_2_z,
                      double pt_3_x, double pt_3_y, double pt_3_z) {
    double a[3], b[3], cross[3], area;
    a[0] = pt_2_x - pt_1_x;
    a[1] = pt_2_y - pt_1_y;
    a[2] = pt_2_z - pt_1_z;
    b[0] = pt_3_x - pt_1_x;
    b[1] = pt_3_y - pt_1_y;
    b[2] = pt_3_z - pt_1_z;
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

    NIMBLE_FUNCTION
    ContactEntity() {}

    NIMBLE_FUNCTION
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
            node_id_1_for_fictitious_node_ = node_ids_for_fictitious_node[0];
            node_id_2_for_fictitious_node_ = node_ids_for_fictitious_node[1];
            node_id_3_for_fictitious_node_ = node_ids_for_fictitious_node[2];
            node_id_4_for_fictitious_node_ = node_ids_for_fictitious_node[3];
          }
          else{
            printf("\n**** Error in ContactEntity constructor, invalid entity type.\n");
          }
          coord_1_x_ = coord[0];
          coord_1_y_ = coord[1];
          coord_1_z_ = coord[2];
          if (entity_type_ == TRIANGLE) {
            coord_2_x_ = coord[3];
            coord_2_y_ = coord[4];
            coord_2_z_ = coord[5];
            coord_3_x_ = coord[6];
            coord_3_y_ = coord[7];
            coord_3_z_ = coord[8];
          }
          SetBoundingBox();
        }

    NIMBLE_FUNCTION
    ~ContactEntity() {}

    template <typename ArgT>
    NIMBLE_INLINE_FUNCTION
    void SetCoordinates(ArgT coord) {
      int n0 = 3*node_id_for_node_1_;
      coord_1_x_ = coord[n0];
      coord_1_y_ = coord[n0+1];
      coord_1_z_ = coord[n0+2];
      if (entity_type_ == TRIANGLE) {
        n0 = 3*node_id_for_node_2_;
        coord_2_x_ = coord[n0];
        coord_2_y_ = coord[n0+1];
        coord_2_z_ = coord[n0+2];
        n0 = 3*node_id_1_for_fictitious_node_;
        int n1 = 3*node_id_2_for_fictitious_node_;
        int n2 = 3*node_id_3_for_fictitious_node_;
        int n3 = 3*node_id_4_for_fictitious_node_;
        coord_3_x_ = (coord[n0]   + coord[n1]   + coord[n2]   + coord[n3]  ) / 4.0;
        coord_3_y_ = (coord[n0+1] + coord[n1+1] + coord[n2+1] + coord[n3+1]) / 4.0;
        coord_3_z_ = (coord[n0+2] + coord[n1+2] + coord[n2+2] + coord[n3+2]) / 4.0;
      }
      SetBoundingBox();
    }

    template <typename ArgT>
    NIMBLE_INLINE_FUNCTION
    void GetForces(ArgT force) {
      int n = 3*node_id_for_node_1_;
      force[n]   += force_1_x_;
      force[n+1] += force_1_y_;
      force[n+2] += force_1_z_;
      if (entity_type_ == TRIANGLE) {
        n = 3*node_id_for_node_2_;
        force[n]   += force_2_x_;
        force[n+1] += force_2_y_;
        force[n+2] += force_2_z_;
        int n = 3*node_id_1_for_fictitious_node_;
        force[n]   += force_3_x_ / 4.0;
        force[n+1] += force_3_y_ / 4.0;
        force[n+2] += force_3_z_ / 4.0;
        n = 3*node_id_2_for_fictitious_node_;
        force[n]   += force_3_x_ / 4.0;
        force[n+1] += force_3_y_ / 4.0;
        force[n+2] += force_3_z_ / 4.0;
        n = 3*node_id_3_for_fictitious_node_;
        force[n]   += force_3_x_ / 4.0;
        force[n+1] += force_3_y_ / 4.0;
        force[n+2] += force_3_z_ / 4.0;
        n = 3*node_id_4_for_fictitious_node_;
        force[n]   += force_3_x_ / 4.0;
        force[n+1] += force_3_y_ / 4.0;
        force[n+2] += force_3_z_ / 4.0;
      }
    }

    // functions required for NimbleSMExtras contact search
    NIMBLE_INLINE_FUNCTION
    double get_x_min() const { return bounding_box_x_min_; }
    NIMBLE_INLINE_FUNCTION
    double get_x_max() const { return bounding_box_x_max_; }
    NIMBLE_INLINE_FUNCTION
    double get_y_min() const { return bounding_box_y_min_; }
    NIMBLE_INLINE_FUNCTION
    double get_y_max() const { return bounding_box_y_max_; }
    NIMBLE_INLINE_FUNCTION
    double get_z_min() const { return bounding_box_z_min_; }
    NIMBLE_INLINE_FUNCTION
    double get_z_max() const { return bounding_box_z_max_; }

    // Functions for bvh contact search
    //int contant_entity_global_id() const { return contant_entity_global_id_; }
    //const double* centroid() const {
    //  // FIX ME
    //  // centroid is now stored as centroid_x_ centroid_y_ and centroid_z_
    //  // using arrays complicates kokkos / cuda mem copies, so I just unrolled everything
    //  return centroid_;
    //}
#ifdef NIMBLE_HAVE_BVH
    bvh::dop_26<double> kdop() const {
      std::vector< std::array< double, 3 > > vertices;
      vertices.push_back( {{ coord_1_x_, coord_1_y_, coord_1_z_ }} );
      if (entity_type_ == TRIANGLE) {
        vertices.push_back( {{ coord_2_x_, coord_2_y_, coord_2_z_ }} );
        vertices.push_back( {{ coord_3_x_, coord_3_y_, coord_3_z_ }} );
      }
      const double inflation_length = 0.15 * char_len_;
      return bvh::dop_26< double >::from_vertices( vertices.begin(), vertices.end(), inflation_length );
    }
#endif

    NIMBLE_INLINE_FUNCTION
    void SetBoundingBox() {

      centroid_x_ = coord_1_x_;
      centroid_y_ = coord_1_y_;
      centroid_z_ = coord_1_z_;
      bounding_box_x_min_ = coord_1_x_;
      bounding_box_x_max_ = coord_1_x_;
      bounding_box_y_min_ = coord_1_y_;
      bounding_box_y_max_ = coord_1_y_;
      bounding_box_z_min_ = coord_1_z_;
      bounding_box_z_max_ = coord_1_z_;

      // todo, try ternary operator here
      if (entity_type_ == TRIANGLE) {
        centroid_x_ += coord_2_x_;
        centroid_y_ += coord_2_y_;
        centroid_z_ += coord_2_z_;
        if (coord_2_x_ < bounding_box_x_min_)
          bounding_box_x_min_ = coord_2_x_;
        if (coord_2_x_ > bounding_box_x_max_)
          bounding_box_x_max_ = coord_2_x_;
        if (coord_2_y_ < bounding_box_y_min_)
          bounding_box_y_min_ = coord_2_y_;
        if (coord_2_y_ > bounding_box_y_max_)
          bounding_box_y_max_ = coord_2_y_;
        if (coord_2_z_ < bounding_box_z_min_)
          bounding_box_z_min_ = coord_2_z_;
        if (coord_2_z_ > bounding_box_z_max_)
          bounding_box_z_max_ = coord_2_z_;
        centroid_x_ += coord_3_x_;
        centroid_y_ += coord_3_y_;
        centroid_z_ += coord_3_z_;
        if (coord_3_x_ < bounding_box_x_min_)
          bounding_box_x_min_ = coord_3_x_;
        if (coord_3_x_ > bounding_box_x_max_)
          bounding_box_x_max_ = coord_3_x_;
        if (coord_3_y_ < bounding_box_y_min_)
          bounding_box_y_min_ = coord_3_y_;
        if (coord_3_y_ > bounding_box_y_max_)
          bounding_box_y_max_ = coord_3_y_;
        if (coord_3_z_ < bounding_box_z_min_)
          bounding_box_z_min_ = coord_3_z_;
        if (coord_3_z_ > bounding_box_z_max_)
          bounding_box_z_max_ = coord_3_z_;
      }

      centroid_x_ /= num_nodes_;
      centroid_y_ /= num_nodes_;
      centroid_z_ /= num_nodes_;

      double inflation_length = 0.15 * char_len_;

      bounding_box_x_min_ -= inflation_length;
      bounding_box_x_max_ += inflation_length;
      bounding_box_y_min_ -= inflation_length;
      bounding_box_y_max_ += inflation_length;
      bounding_box_z_min_ -= inflation_length;
      bounding_box_z_max_ += inflation_length;
    }

    NIMBLE_INLINE_FUNCTION
    void ComputeNodalContactForces(const double * const contact_force,
                                   const double * const closest_point_projection) {

      double N[3] = {0.0, 0.0, 0.0};

      if (entity_type_ == NODE) {
        N[0] = 1.0;
      }
      else if (entity_type_ == TRIANGLE) {
        // Find the natural coordinates of the closest_point_projection
        double area_1 = TriangleArea(closest_point_projection[0], closest_point_projection[1], closest_point_projection[2],
                                     coord_2_x_, coord_2_y_, coord_2_z_,
                                     coord_3_x_, coord_3_y_, coord_3_z_);
        double area_2 = TriangleArea(coord_1_x_, coord_1_y_, coord_1_z_,
                                     closest_point_projection[0], closest_point_projection[1], closest_point_projection[2],
                                     coord_3_x_, coord_3_y_, coord_3_z_);
        double full_area = TriangleArea(coord_1_x_, coord_1_y_, coord_1_z_,
                                        coord_2_x_, coord_2_y_, coord_2_z_,
                                        coord_3_x_, coord_3_y_, coord_3_z_);
        N[0] = area_1/full_area;
        N[1] = area_2/full_area;
        N[2] = 1.0 - N[0] - N[1];
      }

      force_1_x_ = N[0] * contact_force[0];
      force_1_y_ = N[0] * contact_force[1];
      force_1_z_ = N[0] * contact_force[2];
      if (entity_type_ == TRIANGLE) {
        force_2_x_ = N[1] * contact_force[3];
        force_2_y_ = N[1] * contact_force[4];
        force_2_z_ = N[1] * contact_force[5];
        force_3_x_ = N[2] * contact_force[6];
        force_3_y_ = N[2] * contact_force[7];
        force_3_z_ = N[2] * contact_force[8];
      }
    }

    int contant_entity_global_id_ = -1;

    // positions of nodes that define triangular contact patch
    CONTACT_ENTITY_TYPE entity_type_ = NONE;
    int num_nodes_ = 0;
    double coord_1_x_ = 0.0;
    double coord_1_y_ = 0.0;
    double coord_1_z_ = 0.0;
    double coord_2_x_ = 0.0;
    double coord_2_y_ = 0.0;
    double coord_2_z_ = 0.0;
    double coord_3_x_ = 0.0;
    double coord_3_y_ = 0.0;
    double coord_3_z_ = 0.0;
    double force_1_x_ = 0.0;
    double force_1_y_ = 0.0;
    double force_1_z_ = 0.0;
    double force_2_x_ = 0.0;
    double force_2_y_ = 0.0;
    double force_2_z_ = 0.0;
    double force_3_x_ = 0.0;
    double force_3_y_ = 0.0;
    double force_3_z_ = 0.0;

    double char_len_ = 0.0;

    // bounding box for NimbleSMExtras contact search routines
    double bounding_box_x_min_ = -DBL_MAX;
    double bounding_box_x_max_ =  DBL_MAX;
    double bounding_box_y_min_ = -DBL_MAX;
    double bounding_box_y_max_ =  DBL_MAX;
    double bounding_box_z_min_ = -DBL_MAX;
    double bounding_box_z_max_ =  DBL_MAX;

    // centroid for bvh search routines
    double centroid_x_ = 0.0;
    double centroid_y_ = 0.0;
    double centroid_z_ = 0.0;

    // map for moving displacement/forces to/from the contact manager data structures
    int node_id_for_node_1_ = -1;
    int node_id_for_node_2_ = -1;
    // fictitious node maps to multiple real nodes
    int node_id_1_for_fictitious_node_ = -1;
    int node_id_2_for_fictitious_node_ = -1;
    int node_id_3_for_fictitious_node_ = -1;
    int node_id_4_for_fictitious_node_ = -1;
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

      int num_nodes_in_contact_manager = node_ids_d_.extent(0);

      Kokkos::parallel_for("ContactManager::ApplyDisplacements set coord_d_ vector",
                           num_nodes_in_contact_manager,
                           KOKKOS_LAMBDA(const int i) {
        //int node_id = node_ids_d_(i);
        //printf("\nDEBUGGING node_id %d", node_id);
        //coord_d_(3*i)   = model_coord_d_(3*i)   + displacement_d(node_id, 0);
        //coord_d_(3*i+1) = model_coord_d_(3*i+1) + displacement_d(node_id, 1);
        //coord_d_(3*i+2) = model_coord_d_(3*i+2) + displacement_d(node_id, 2);
      });
      /*
      Kokkos::parallel_for("ContactManager::ApplyDisplacements set contact node entity displacements",
                         num_nodes_in_contact_manager,
                         KOKKOS_LAMBDA(const int i_node) {
        contact_nodes_d_(i_node).SetCoordinates(coord_d_);
      });

      Kokkos::parallel_for("ContactManager::ApplyDisplacements set contact face entity displacements",
                         num_nodes_in_contact_manager,
                         KOKKOS_LAMBDA(const int i_face) {
        contact_faces_d_(i_face).SetCoordinates(coord_d_);
      });
      */
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
      /*
      Kokkos::parallel_for("ContactManager::GetForces",
                         num_nodes_in_contact_manager,
                         KOKKOS_LAMBDA(const int i) {
        int node_id = node_ids_d_(i);
        contact_force_d(node_id, 0) = force_d_(3*i);
        contact_force_d(node_id, 1) = force_d_(3*i+1);
        contact_force_d(node_id, 2) = force_d_(3*i+2);
      });
      */
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
