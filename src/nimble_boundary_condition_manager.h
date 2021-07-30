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

#ifndef NIMBLE_BOUNDARY_CONDITION_MANAGER_H
#define NIMBLE_BOUNDARY_CONDITION_MANAGER_H

#include "nimble_boundary_condition.h"
#include "nimble_view.h"

#ifdef NIMBLE_HAVE_DARMA
#include "darma.h"
#else
#include <map>
#include <vector>
#endif

namespace nimble {

class BoundaryConditionManager
{
 public:
  enum Time_Integration_Scheme
  {
    EXPLICIT    = 0,
    QUASISTATIC = 1,
    UNDEFINED   = 2
  };

  BoundaryConditionManager() {}

  void
  Initialize(
      std::map<int, std::string> const&      node_set_names,
      std::map<int, std::vector<int>> const& node_sets,
      std::map<int, std::string> const&      side_set_names,
      std::map<int, std::vector<int>> const& side_sets,
      std::vector<std::string> const&        bc_strings,
      int                                    dim,
      std::string                            time_integration_scheme);

  virtual ~BoundaryConditionManager() = default;

#ifdef NIMBLE_HAVE_DARMA
  template <typename ArchiveType>
  void
  serialize(ArchiveType& ar)
  {
    ar | node_set_names_ | node_sets_ | side_set_names_ | side_sets_ | boundary_conditions_ | dim_ |
        time_integration_scheme_;
  }
#endif

  template <typename ViewT>
  void
  ApplyInitialConditions(const ViewT reference_coordinates, ViewT velocity, std::vector<ViewT>& offnom_velocities)
  {
    for (auto& bc : boundary_conditions_) {
      if (bc.bc_type_ == BoundaryCondition::INITIAL_VELOCITY) {
        int                     node_set_id    = bc.node_set_id_;
        int                     coordinate     = bc.coordinate_;
        bool                    has_expression = bc.has_expression_;
        std::vector<int> const& node_set       = node_sets_[node_set_id];
        if (!has_expression) {
          double magnitude = bc.magnitude_;
          for (int n : node_set) {
            velocity(n, coordinate) = magnitude;
            for (int nuq = 0; nuq < offnom_velocities.size(); nuq++) {
              offnom_velocities[nuq](n, coordinate) = magnitude;
            }
          }
        } else {
          for (int n : node_set) {
            bc.expression_.x = reference_coordinates(n, 0);
            bc.expression_.y = reference_coordinates(n, 1);
            bc.expression_.z = 0.0;
            if (dim_ == 3) { bc.expression_.z = reference_coordinates(n, 2); }
            bc.expression_.t        = 0.0;
            velocity(n, coordinate) = bc.expression_.eval();
            for (int nuq = 0; nuq < offnom_velocities.size(); nuq++) {
              offnom_velocities[nuq](n, coordinate) = bc.expression_.eval();
            }
          }
        }
      }
    }
  }

  template <typename ViewT>
  void
  ApplyInitialConditions(const ViewT reference_coordinates, ViewT velocity)
  {
    std::vector<ViewT> empty;
    ApplyInitialConditions<ViewT>(reference_coordinates, velocity, empty);
  }

  template <typename ViewT>
  void
  ApplyKinematicBC(
      double              time_current,
      double              time_previous,
      const ViewT         reference_coordinates,
      ViewT               displacement,
      ViewT               velocity,
      std::vector<ViewT>& offnom_velocities)
  {
    double delta_t = time_current - time_previous;

    for (auto& bc : boundary_conditions_) {
      int                     node_set_id    = bc.node_set_id_;
      int                     coordinate     = bc.coordinate_;
      bool                    has_expression = bc.has_expression_;
      std::vector<int> const& node_set       = node_sets_[node_set_id];

      if (bc.bc_type_ == BoundaryCondition::PRESCRIBED_VELOCITY) {
        if (!has_expression) {
          double velocity_magnitude = bc.magnitude_;
          for (int n : node_set) {
            velocity(n, coordinate) = velocity_magnitude;
            for (int nuq = 0; nuq < offnom_velocities.size(); nuq++) {
              offnom_velocities[nuq](n, coordinate) = velocity_magnitude;
            }
            if (time_integration_scheme_ == QUASISTATIC) {
              displacement(n, coordinate) += velocity_magnitude * delta_t;
            }
          }
        } else {
          for (int n : node_set) {
            bc.expression_.x = reference_coordinates(n, 0);
            bc.expression_.y = reference_coordinates(n, 1);
            bc.expression_.z = 0.0;
            if (dim_ == 3) { bc.expression_.z = reference_coordinates(n, 2); }
            bc.expression_.t          = time_current;
            double velocity_magnitude = bc.expression_.eval();
            velocity(n, coordinate)   = velocity_magnitude;
            for (int nuq = 0; nuq < offnom_velocities.size(); nuq++) {
              offnom_velocities[nuq](n, coordinate) = velocity_magnitude;
            }
            if (time_integration_scheme_ == QUASISTATIC) {
              displacement(n, coordinate) += velocity_magnitude * delta_t;
            }
          }
        }
      } else if (bc.bc_type_ == BoundaryCondition::PRESCRIBED_DISPLACEMENT && delta_t > 0.0) {
        if (!has_expression) {
          double displacement_magnitude = bc.magnitude_;
          for (int n : node_set) {
            velocity(n, coordinate) = (displacement_magnitude - displacement(n, coordinate)) / delta_t;
            if (time_integration_scheme_ == QUASISTATIC) { displacement(n, coordinate) = displacement_magnitude; }
          }
        } else {
          for (int n : node_set) {
            bc.expression_.x = reference_coordinates(n, 0);
            bc.expression_.y = reference_coordinates(n, 1);
            bc.expression_.z = 0.0;
            if (dim_ == 3) { bc.expression_.z = reference_coordinates(n, 2); }
            bc.expression_.t              = time_current;
            double displacement_magnitude = bc.expression_.eval();
            velocity(n, coordinate)       = (displacement_magnitude - displacement(n, coordinate)) / delta_t;
            if (time_integration_scheme_ == QUASISTATIC) { displacement(n, coordinate) = displacement_magnitude; }
          }
        }
      }
    }
  }

  template <typename ViewT>
  void
  ApplyKinematicBC(
      double      time_current,
      double      time_previous,
      const ViewT reference_coordinates,
      ViewT       displacement,
      ViewT       velocity)
  {
    std::vector<ViewT> empty;
    ApplyKinematicBC(time_current, time_previous, reference_coordinates, displacement, velocity, empty);
  }

  template <typename MatT>
  void
  ModifyTangentStiffnessMatrixForKinematicBC(
      int        num_unknowns,
      const int* global_node_ids,
      double     diagonal_entry,
      MatT&      tangent_stiffness) const;

  void
  ModifyRHSForKinematicBC(const int* global_node_ids, double* rhs) const;

  template <typename ViewT>
  void
  ApplyTractionBC(double time_current, const ViewT reference_coordinates, ViewT internal_force)
  {
    for (auto& bc : boundary_conditions_) {
      int                     side_set_id    = bc.side_set_id_;
      int                     coordinate     = bc.coordinate_;
      bool                    has_expression = bc.has_expression_;
      std::vector<int> const& side_set       = side_sets_[side_set_id];
      if (bc.bc_type_ == BoundaryCondition::PRESCRIBED_TRACTION) {
        if (has_expression == true) {
          for (int ss : side_set) {}
        } else {
          for (int ss : side_set) {}
        }
      }
    }
  }

 private:
  std::map<int, std::string>      node_set_names_;
  std::map<int, std::vector<int>> node_sets_;
  std::map<int, std::string>      side_set_names_;
  std::map<int, std::vector<int>> side_sets_;
  std::vector<BoundaryCondition>  boundary_conditions_;
  int                             dim_{0};
  Time_Integration_Scheme         time_integration_scheme_{UNDEFINED};
};

}  // namespace nimble

#endif  // NIMBLE_BOUNDARY_CONDITION_MANAGER_H
