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

  NIMBLE_INLINE_FUNCTION
  double PointEdgeClosestPointFindT(double const p1[],
                                    double const p2[],
                                    double const p[]) {
    return ((p[0] - p1[0] )*(p2[0] - p1[0]) + (p[1] - p1[1] )*(p2[1] - p1[1]) + (p[2] - p1[2] )*(p2[2] - p1[2]))
      / ((p2[0] - p1[0])*(p2[0] - p1[0]) + (p2[1] - p1[1])*(p2[1] - p1[1]) + (p2[2] - p1[2])*(p2[2] - p1[2]));
  }

  NIMBLE_INLINE_FUNCTION
  double PointEdgeClosestPointFindDistanceSquared(double const p1[],
                                                  double const p2[],
                                                  double const p[],
                                                  double t) {
    double temp1 = p1[0] + (p2[0] - p1[0])*t - p[0];
    double temp2 = p1[1] + (p2[1] - p1[1])*t - p[1];
    double temp3 = p1[2] + (p2[2] - p1[2])*t - p[2];
    return (temp1*temp1 + temp2*temp2 + temp3*temp3);
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

    struct vertex {
      vertex() {
        coords_[0] = coords_[1] = coords_[2] = 0.0;
      }
      double& operator[](int i) {
        return coords_[i];
      }
      double operator[](int i) const {
        return coords_[i];
      }
      double coords_[3];
    };

    NIMBLE_FUNCTION
    ContactEntity() {}

    NIMBLE_FUNCTION
    ContactEntity(CONTACT_ENTITY_TYPE entity_type,
                  int contact_entity_global_id,
                  double const coord[],
                  double characteristic_length,
                  int node_id_for_node_1,
                  int node_id_for_node_2 = 0,
                  int node_ids_for_fictitious_node[4] = 0)
      : entity_type_(entity_type),
        contact_entity_global_id_(contact_entity_global_id),
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
    int contact_entity_global_id() const { return contact_entity_global_id_; }

    vertex centroid() const {
      return centroid_;
    }

#ifdef NIMBLE_HAVE_BVH
    bvh::dop_26<double> kdop() const {
      const double inflation_length = 0.15 * char_len_;
      if (entity_type_ == NODE) {
        vertex v;
        v[0] = coord_1_x_;
        v[1] = coord_1_y_;
        v[2] = coord_1_z_;
        return bvh::dop_26< double >::from_vertices( &v, &v + 1, inflation_length );
      }
      // entity_type_ == TRIANGLE
      vertex v[3];
      v[0][0] = coord_1_x_;
      v[0][1] = coord_1_y_;
      v[0][2] = coord_1_z_;
      v[1][0] = coord_2_x_;
      v[1][1] = coord_2_y_;
      v[1][2] = coord_2_z_;
      v[2][0] = coord_3_x_;
      v[2][1] = coord_3_y_;
      v[2][2] = coord_3_z_;
      return bvh::dop_26< double >::from_vertices( v, v + 3, inflation_length );
    }
#endif

    NIMBLE_INLINE_FUNCTION
    void SetBoundingBox() {

      centroid_[0] = coord_1_x_;
      centroid_[1] = coord_1_y_;
      centroid_[2] = coord_1_z_;
      bounding_box_x_min_ = coord_1_x_;
      bounding_box_x_max_ = coord_1_x_;
      bounding_box_y_min_ = coord_1_y_;
      bounding_box_y_max_ = coord_1_y_;
      bounding_box_z_min_ = coord_1_z_;
      bounding_box_z_max_ = coord_1_z_;

      // todo, try ternary operator here
      if (entity_type_ == TRIANGLE) {
        centroid_[0] += coord_2_x_;
        centroid_[1] += coord_2_y_;
        centroid_[2] += coord_2_z_;
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
        centroid_[0] += coord_3_x_;
        centroid_[1] += coord_3_y_;
        centroid_[2] += coord_3_z_;
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

    int contact_entity_global_id_ = -1;

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
    vertex centroid_;

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

    typedef enum ProjectionType {
      UNKNOWN=0,
      NODE_OR_EDGE=1,
      FACE=2
    } PROJECTION_TYPE;

    ContactManager() : penalty_parameter_(0.0)
#ifdef NIMBLE_HAVE_EXTRAS
      , contact_nodes_search_tree_("contact nodes search tree")
      , contact_faces_search_tree_("contact faces search tree")
#endif
      {}

    virtual ~ContactManager() {}

    bool ContactEnabled() { return contact_enabled_; }

    std::vector< std::vector<int> > SkinBlocks(GenesisMesh const & mesh,
                                               std::vector<int> const & partition_boundary_node_ids,
                                               std::vector<int> const & block_ids);

    void SetPenaltyParameter(double penalty_parameter) {
      penalty_parameter_ = penalty_parameter;
    }

    void CreateContactEntities(GenesisMesh const & mesh,
                               std::vector<int> const & partition_boundary_node_ids,
                               std::vector<int> const & master_block_ids,
                               std::vector<int> const & slave_block_ids);

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
      int num_contact_node_entities = contact_nodes_d_.extent(0);
      int num_contact_face_entities = contact_faces_d_.extent(0);

      // circumvent lambda *this glitch
      nimble_kokkos::DeviceIntegerArrayView node_ids = node_ids_d_;
      nimble_kokkos::DeviceScalarNodeView model_coord = model_coord_d_;
      nimble_kokkos::DeviceScalarNodeView coord = coord_d_;
      DeviceContactEntityArrayView contact_nodes = contact_nodes_d_;
      DeviceContactEntityArrayView contact_faces = contact_faces_d_;

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

    void ClosestPointProjection(std::vector<ContactEntity> const & nodes,
                                std::vector<ContactEntity> const & triangles,
                                std::vector<ContactEntity::vertex>& closest_points,
                                std::vector<PROJECTION_TYPE>& projection_types);

    void ComputeContactForce(int step, bool debug_output);

    void WriteContactEntitiesToVTKFile(int step);

    void WriteContactEntitiesToVTKFile(const std::vector<ContactEntity> &faces,
                                       const std::vector<ContactEntity> &nodes,
                                       const std::string &prefix,
                                       int step);

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

    using DeviceContactEntityArrayView = Kokkos::View< ContactEntity*, nimble_kokkos::kokkos_layout, nimble_kokkos::kokkos_device >;
    DeviceContactEntityArrayView contact_faces_d_ = DeviceContactEntityArrayView("contact_faces_d", 1);
    DeviceContactEntityArrayView contact_nodes_d_ = DeviceContactEntityArrayView("contact_nodes_d", 1);

    // TODO remove this once enforcement is on device
    using HostContactEntityArrayView = Kokkos::View< ContactEntity*, nimble_kokkos::kokkos_layout, nimble_kokkos::kokkos_host >;
    HostContactEntityArrayView contact_faces_h_ = HostContactEntityArrayView("contact_faces_h", 1);
    HostContactEntityArrayView contact_nodes_h_ = HostContactEntityArrayView("contact_nodes_h", 1);
#endif

#ifdef NIMBLE_HAVE_EXTRAS
    stk::search::mas_aabb_tree<double, nimble_kokkos::kokkos_device_execution_space> contact_nodes_search_tree_;
    stk::search::mas_aabb_tree<double, nimble_kokkos::kokkos_device_execution_space> contact_faces_search_tree_;
#endif
  };

} // namespace nimble


#ifdef NIMBLE_HAVE_BVH
namespace bvh
{
  /**
   * Adapter for bvh to get contact entity collision info
   */
  template<>
  struct get_entity_info< nimble::ContactEntity >
  {
    static decltype(auto) get_kdop( const nimble::ContactEntity &_entity )
    {
      return _entity.kdop();
    }
    
    static decltype(auto) get_global_id( const nimble::ContactEntity &_entity )
    {
      return _entity.contact_entity_global_id();
    }
    
    static decltype(auto) get_centroid( const nimble::ContactEntity &_entity )
    {
      return _entity.centroid();
    }
  };
}
#endif  // NIMBLE_HAVE_BVH

#endif // NIMBLE_MATERIAL_H
