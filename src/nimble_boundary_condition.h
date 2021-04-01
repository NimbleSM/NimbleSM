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

#ifndef NIMBLE_BOUNDARY_CONDITION_H
#define NIMBLE_BOUNDARY_CONDITION_H

#include "nimble_expression_parser.h"

#ifdef NIMBLE_HAVE_DARMA
#include "darma.h"
#else
#include <map>
#include <string>
#include <vector>
#endif

namespace nimble {

class BoundaryCondition
{
 public:
  enum Boundary_Condition_Type
  {
    INITIAL_VELOCITY        = 0,
    PRESCRIBED_VELOCITY     = 1,
    PRESCRIBED_DISPLACEMENT = 2,
    PERIODIC_RVE            = 3,
    RVE_FIXED_DISPLACEMENT  = 4,
    UNDEFINED               = 5
  };

  BoundaryCondition()
      : dim_(0),
        node_set_name_("unknown"),
        node_set_id_(-1),
        coordinate_(-1),
        magnitude_(0.0),
        bc_type_(UNDEFINED),
        has_expression_(false)
  {
  }

  virtual ~BoundaryCondition() {}

#ifdef NIMBLE_HAVE_DARMA
  template <typename ArchiveType>
  void
  serialize(ArchiveType& ar)
  {
    ar | dim_ | node_set_name_ | node_set_id_ | coordinate_ | magnitude_ |
        bc_type_ | has_expression_ | expression_string_ |
        rve_macroscale_deformation_gradient_strings_;

    if (ar.is_unpacking()) {
      if (has_expression_) {
        expression_ =
            ExpressionParsing::BoundaryConditionFunctor(expression_string_);
      }
    }
  }
#endif

  bool
  Initialize(
      int                               dim,
      std::string                       bc_string,
      std::map<int, std::string> const& node_set_names);

  void
  GetRVEMacroscaleDeformationGradient(
      double  time,
      double* deformation_gradient,
      double  x,
      double  y,
      double  z);

  int                                         dim_;
  std::string                                 node_set_name_;
  int                                         node_set_id_;
  int                                         coordinate_;
  double                                      magnitude_;
  Boundary_Condition_Type                     bc_type_;
  bool                                        has_expression_;
  std::string                                 expression_string_;
  ExpressionParsing::BoundaryConditionFunctor expression_;
  std::vector<std::string> rve_macroscale_deformation_gradient_strings_;
};

}  // namespace nimble

#endif  // NIMBLE_BOUNDARY_CONDITION_H
