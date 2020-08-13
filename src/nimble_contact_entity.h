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

#ifndef NIMBLE_CONTACTENTITY_H
#define NIMBLE_CONTACTENTITY_H

#include <vector>
#include <cstring>
#include <string>
#include <cfloat>
#include <math.h>
#include "nimble_kokkos_defs.h"
#include "nimble_utils.h"

#ifdef NIMBLE_HAVE_BVH
  #include <bvh/kdop.hpp>
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
    CrossProduct(b, a, cross);
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
      NIMBLE_INLINE_FUNCTION
      vertex() {
        coords_[0] = coords_[1] = coords_[2] = 0.0;
      }
      NIMBLE_INLINE_FUNCTION
      double& operator[](int i) {
        return coords_[i];
      }
      NIMBLE_INLINE_FUNCTION
      double operator[](int i) const {
        return coords_[i];
      }
      double coords_[3];
    };

    NIMBLE_INLINE_FUNCTION
    ContactEntity() = default;

    NIMBLE_INLINE_FUNCTION
    ContactEntity(CONTACT_ENTITY_TYPE entity_type,
                  int contact_entity_global_id,
                  double const coord[],
                  double characteristic_length,
                  int node_id_for_node_1,
                  int node_id_for_node_2 = 0,
                  const int node_ids_for_fictitious_node[4] = 0)
      : entity_type_(entity_type),
        char_len_(characteristic_length),
        contact_entity_global_id_(contact_entity_global_id){

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

    NIMBLE_INLINE_FUNCTION
    ~ContactEntity() = default;

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
    void ScatterForceToContactManagerForceVector(ArgT force) {
      int n = 3*node_id_for_node_1_;
      force[n]   += force_1_x_;
      force[n+1] += force_1_y_;
      force[n+2] += force_1_z_;
      if (entity_type_ == TRIANGLE) {
        n = 3*node_id_for_node_2_;
        force[n]   += force_2_x_;
        force[n+1] += force_2_y_;
        force[n+2] += force_2_z_;
        n = 3*node_id_1_for_fictitious_node_;
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
    NIMBLE_INLINE_FUNCTION
    int contact_entity_global_id() const { return contact_entity_global_id_; }

    NIMBLE_INLINE_FUNCTION
    vertex centroid() const {
      return centroid_;
    }

    NIMBLE_INLINE_FUNCTION
    bool contact_status() const noexcept {
      return contact_status_;
    }

    NIMBLE_INLINE_FUNCTION
    void set_contact_status( bool status ) noexcept {
      contact_status_ = status;
    }

#ifdef NIMBLE_HAVE_BVH
    bvh::dop_26d kdop_;

    void RecomputeKdop()
    {
      const double inflation_length = inflation_factor * char_len_;
      if (entity_type_ == NODE) {
        vertex v;
        v[0] = coord_1_x_;
        v[1] = coord_1_y_;
        v[2] = coord_1_z_;
        kdop_ = bvh::dop_26d::from_vertices( &v, &v + 1, inflation_length );
      } else {
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
        kdop_ = bvh::dop_26d::from_vertices(v, v + 3, inflation_length);
      }
    }

    bvh::dop_26<double> Kdop() const {
      return kdop_;
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

      double inflation_length = inflation_factor * char_len_;

      bounding_box_x_min_ -= inflation_length;
      bounding_box_x_max_ += inflation_length;
      bounding_box_y_min_ -= inflation_length;
      bounding_box_y_max_ += inflation_length;
      bounding_box_z_min_ -= inflation_length;
      bounding_box_z_max_ += inflation_length;
    }

    NIMBLE_INLINE_FUNCTION
    void SetNodalContactForces(const double * const contact_force,
                               const double * const N = NULL) {
      if (entity_type_ == NODE) {
        force_1_x_ = -contact_force[0];
        force_1_y_ = -contact_force[1];
        force_1_z_ = -contact_force[2];
      }
      else if (entity_type_ == TRIANGLE) {
        force_1_x_ = N[0] * contact_force[0];
        force_1_y_ = N[0] * contact_force[1];
        force_1_z_ = N[0] * contact_force[2];
        force_2_x_ = N[1] * contact_force[0];
        force_2_y_ = N[1] * contact_force[1];
        force_2_z_ = N[1] * contact_force[2];
        force_3_x_ = N[2] * contact_force[0];
        force_3_y_ = N[2] * contact_force[1];
        force_3_z_ = N[2] * contact_force[2];
      }
    }

// DEPRECATED
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
        force_2_x_ = N[1] * contact_force[0];
        force_2_y_ = N[1] * contact_force[1];
        force_2_z_ = N[1] * contact_force[2];
        force_3_x_ = N[2] * contact_force[0];
        force_3_y_ = N[2] * contact_force[1];
        force_3_z_ = N[2] * contact_force[2];
      }
    }

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

    /// \brief Factor multiplying characteristic length to inflate bounding box
    ///
    /// \note Empirical default value 0.15
    double inflation_factor = 0.15;

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

    // store a unique global id
    // NODE:
    //   the contact_entity_global_id_ is the exodus global node id in the parent FEM mesh
    // FACE:
    //   the contact_entity_global_id_ is a bit-wise combination of the global element id, the face ordinal, and the triangle ordinal
    //   first 2 bits are the triangle ordinal (range is 1-4)
    //   next 3 bits are the face ordinal (range is 1-6)
    //   remaining bits are the genesis element id from the parent FEM mesh (e.g., the global id of the hex from which the face was extracted)
    int contact_entity_global_id_ = -1;

    bool contact_status_ = false;
  };

  template <typename ArgT>
  void SerializeContactFaces(int num_entities,
                             ArgT contact_entities,
                             std::vector<char>& buffer) {
    constexpr size_t size_int = sizeof(int);
    constexpr size_t size_double = sizeof(double);
    size_t entity_size = 10*sizeof(double) + 8*sizeof(int);
    buffer.resize(entity_size*num_entities);
    char *scan = &buffer[0];
    for (int i=0 ; i<num_entities ; i++) {
      memcpy(scan, &contact_entities[i].coord_1_x_, size_double); scan += size_double;
      memcpy(scan, &contact_entities[i].coord_1_y_, size_double); scan += size_double;
      memcpy(scan, &contact_entities[i].coord_1_z_, size_double); scan += size_double;
      memcpy(scan, &contact_entities[i].coord_2_x_, size_double); scan += size_double;
      memcpy(scan, &contact_entities[i].coord_2_y_, size_double); scan += size_double;
      memcpy(scan, &contact_entities[i].coord_2_z_, size_double); scan += size_double;
      memcpy(scan, &contact_entities[i].coord_3_x_, size_double); scan += size_double;
      memcpy(scan, &contact_entities[i].coord_3_y_, size_double); scan += size_double;
      memcpy(scan, &contact_entities[i].coord_3_z_, size_double); scan += size_double;
      memcpy(scan, &contact_entities[i].char_len_, size_double); scan += size_double;
      memcpy(scan, &contact_entities[i].contact_entity_global_id_, size_int); scan += size_int;
      memcpy(scan, &contact_entities[i].node_id_for_node_1_, size_int); scan += size_int;
      memcpy(scan, &contact_entities[i].node_id_for_node_2_, size_int); scan += size_int;
      memcpy(scan, &contact_entities[i].node_id_1_for_fictitious_node_, size_int); scan += size_int;
      memcpy(scan, &contact_entities[i].node_id_2_for_fictitious_node_, size_int); scan += size_int;
      memcpy(scan, &contact_entities[i].node_id_3_for_fictitious_node_, size_int); scan += size_int;
      memcpy(scan, &contact_entities[i].node_id_4_for_fictitious_node_, size_int); scan += size_int;
    }
  }

  template <typename ArgT>
  void UnserializeContactFaces(int num_entities,
                               ArgT contact_entities,
                               std::vector<char>& buffer) {
    constexpr size_t size_int = sizeof(int);
    constexpr size_t size_double = sizeof(double);
    ContactEntity entity;
    entity.entity_type_ = ContactEntity::TRIANGLE;
    entity.num_nodes_ = 3;
    entity.force_1_x_ = 0.0;
    entity.force_1_y_ = 0.0;
    entity.force_1_z_ = 0.0;
    entity.force_2_x_ = 0.0;
    entity.force_2_y_ = 0.0;
    entity.force_2_z_ = 0.0;
    entity.force_3_x_ = 0.0;
    entity.force_3_y_ = 0.0;
    entity.force_3_z_ = 0.0;
    char *scan = &buffer[0];
    for (int i=0 ; i<num_entities ; i++) {
      memcpy(&entity.coord_1_x_, scan, size_double); scan += size_double;
      memcpy(&entity.coord_1_y_, scan, size_double); scan += size_double;
      memcpy(&entity.coord_1_z_, scan, size_double); scan += size_double;
      memcpy(&entity.coord_2_x_, scan, size_double); scan += size_double;
      memcpy(&entity.coord_2_y_, scan, size_double); scan += size_double;
      memcpy(&entity.coord_2_z_, scan, size_double); scan += size_double;
      memcpy(&entity.coord_3_x_, scan, size_double); scan += size_double;
      memcpy(&entity.coord_3_y_, scan, size_double); scan += size_double;
      memcpy(&entity.coord_3_z_, scan, size_double); scan += size_double;
      memcpy(&entity.char_len_, scan, size_double); scan += size_double;
      memcpy(&entity.contact_entity_global_id_, scan, size_int); scan += size_int;
      memcpy(&entity.node_id_for_node_1_, scan, size_int); scan += size_int;
      memcpy(&entity.node_id_for_node_2_, scan, size_int); scan += size_int;
      memcpy(&entity.node_id_1_for_fictitious_node_, scan, size_int); scan += size_int;
      memcpy(&entity.node_id_2_for_fictitious_node_, scan, size_int); scan += size_int;
      memcpy(&entity.node_id_3_for_fictitious_node_, scan, size_int); scan += size_int;
      memcpy(&entity.node_id_4_for_fictitious_node_, scan, size_int); scan += size_int;
      entity.SetBoundingBox();
      contact_entities[i] = entity;
    }
  }

#ifdef NIMBLE_HAVE_BVH
  // Free functions for accessing entity info for bvh
  bvh::dop_26<double> get_entity_kdop( const ContactEntity &_entity );
  std::size_t get_entity_global_id( const ContactEntity &_entity );
  bvh::m::vec3d get_entity_centroid( const ContactEntity &_entity );
#endif

} // namespace nimble

#endif // NIMBLE_CONTACTENTITY_H
