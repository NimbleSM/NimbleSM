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
#include "nimble_exodus_output.h"
#include "nimble_kokkos_defs.h"
#include "nimble.mpi.utils.h"

#ifdef NIMBLE_HAVE_EXTRAS
  #include "nimble_extras_contact_includes.h"
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
    cross[1] = b[2]*a[0] - b[0]*a[2];
    cross[2] = b[0]*a[1] - b[1]*a[0];
    area = 0.5 * sqrt(cross[0]*cross[0] + cross[1]*cross[1] + cross[2]*cross[2]);
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
    ContactEntity() {}

    NIMBLE_INLINE_FUNCTION
    ContactEntity(CONTACT_ENTITY_TYPE entity_type,
                  int contact_entity_global_id,
                  double const coord[],
                  double characteristic_length,
                  int node_id_for_node_1,
                  int node_id_for_node_2 = 0,
                  int node_ids_for_fictitious_node[4] = 0,
                  int face_id = -1)
      : entity_type_(entity_type),
        contact_entity_global_id_(contact_entity_global_id),
        char_len_(characteristic_length),
        face_id_(face_id){

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

#ifdef NIMBLE_HAVE_BVH
    bvh::dop_26d kdop_;
    
    void RecomputeKdop()
    {
      const double inflation_length = 0.15 * char_len_;
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
        force_2_x_ = N[1] * contact_force[0];
        force_2_y_ = N[1] * contact_force[1];
        force_2_z_ = N[1] * contact_force[2];
        force_3_x_ = N[2] * contact_force[0];
        force_3_y_ = N[2] * contact_force[1];
        force_3_z_ = N[2] * contact_force[2];
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

    // for faces, store a unique identifier that can be used as a comparision tiebreaker
    // store the global element id, face ordinal, and triangle ordinal as a single int
    // first 2 bits are the triangle ordinal (range is 1-4)
    // next 3 bits are the face ordinal (range is 1-6)
    // remaining bits are the genesis element id
    int face_id_;
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
      memcpy(scan, &contact_entities[i].face_id_, size_int); scan += size_int;
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
      memcpy(&entity.face_id_, scan, size_int); scan += size_int;
      entity.SetBoundingBox();
      contact_entities[i] = entity;
    }
  }
#ifdef NIMBLE_HAVE_BVH
  // Free functions for accessing entity info for bvh
  bvh::dop_26<double> get_entity_kdop( const ContactEntity &_entity );
  std::size_t get_entity_global_id( const ContactEntity &_entity );
  tim::vec3d get_entity_centroid( const ContactEntity &_entity );
#endif

  class ContactManager {

  public:

    typedef enum ProjectionType {
      UNKNOWN=0,
      NODE_OR_EDGE=1,
      FACE=2
    } PROJECTION_TYPE;
    
    explicit ContactManager( std::size_t dicing_factor = 1);

    virtual ~ContactManager() {}

    bool ContactEnabled() { return contact_enabled_; }

    void SkinBlocks(GenesisMesh const & mesh,
                    std::vector<int> const & block_ids,
                    std::vector< std::vector<int> > & skin_faces,
                    std::vector<int> & face_ids);

    void SetPenaltyParameter(double penalty_parameter) {
      penalty_parameter_ = penalty_parameter;
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
                                    std::vector<int> const & master_skin_face_ids,
                                    std::vector<int> const & slave_node_ids,
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
    
    void ComputeContactForce(int step, bool debug_output);

    void ClosestPointProjection(std::vector<ContactEntity> const & nodes,
                                std::vector<ContactEntity> const & triangles,
                                std::vector<ContactEntity::vertex>& closest_points,
                                std::vector<PROJECTION_TYPE>& projection_types);
    
#if defined(NIMBLE_HAVE_MPI) && defined(NIMBLE_HAVE_BVH)
    using patch_collection = bvh::vt::collection< bvh::patch< ContactEntity >, bvh::vt::index_1d >;
    
    void ComputeParallelContactForce(int step, bool is_output_step, bool visualize = false);
#endif

    void InitializeContactVisualization(std::string const & contact_visualization_exodus_file_name);

    void ContactVisualizationWriteStep(double time_current);

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

    double contact_visualization_model_coord_bounding_box_[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    nimble::GenesisMesh genesis_mesh_for_contact_visualization_;
    nimble::ExodusOutput exodus_output_for_contact_visualization_;

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
    
    std::size_t dicing_factor_; ///< Patch dicing factor for overdecomposition

#if defined(NIMBLE_HAVE_MPI) && defined(NIMBLE_HAVE_BVH)
    
    bvh::vt::collision_world<bvh::patch<ContactEntity>, bvh::bvh_tree_26d>  collision_world_;
    
    patch_collection face_patch_collection_;
    patch_collection node_patch_collection_;
#endif

public:

#ifdef NIMBLE_HAVE_EXTRAS
    void compute_and_scatter_contact_force(DeviceContactEntityArrayView contact_nodes_d,
                                           DeviceContactEntityArrayView contact_faces_d,
                                           stk::search::CollisionList<nimble_kokkos::kokkos_device_execution_space> collision_list,
                                           nimble_kokkos::DeviceScalarNodeView contact_manager_force_d);

    void compute_and_scatter_contact_force_OLD(stk::search::CollisionList<nimble_kokkos::kokkos_device_execution_space> collision_list,
                                               DeviceContactEntityArrayView contact_nodes_d,
                                               DeviceContactEntityArrayView contact_faces_d,
                                               int num_contact_nodes,
                                               nimble_kokkos::DeviceScalarNodeView contact_manager_force_d);

    void load_contact_points_and_tris(const int num_collisions,
                                      stk::search::CollisionList<nimble_kokkos::kokkos_device_execution_space> collision_list,
                                      DeviceContactEntityArrayView contact_nodes_d,
                                      DeviceContactEntityArrayView contact_faces_d,
                                      gtk::PointsView<nimble_kokkos::kokkos_device_execution_space> points,
                                      gtk::TrianglesView<nimble_kokkos::kokkos_device_execution_space> triangles);

    KOKKOS_FUNCTION void scatter_contact_forces(double gap,
                                                double scale,
                                                double tri_normal[3],
                                                DeviceContactEntityArrayView contact_faces_d,
                                                int contact_face_index,
                                                const double* closest_pt,
                                                DeviceContactEntityArrayView contact_nodes_d,
                                                int contact_node_index,
                                                nimble_kokkos::DeviceScalarNodeView contact_manager_force_d);
#endif
};
} // namespace nimble

#endif // NIMBLE_MATERIAL_H
