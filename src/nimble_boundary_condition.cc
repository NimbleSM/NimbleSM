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

#include "nimble_boundary_condition.h"

#include <algorithm>
#include <sstream>
#include <stdexcept>

namespace nimble {

namespace {
int
find_name_id(std::map<int, std::string> const& names_map, std::string const& name)
{
  for (auto const& kv : names_map) {
    auto const id       = kv.first;
    auto const map_name = kv.second;
    if (map_name == name) return id;
  }
  return -1;
}
}  // namespace

bool
BoundaryCondition::Initialize(
    int                               dim,
    std::string                       bc_string,
    std::map<int, std::string> const& node_set_names,
    std::map<int, std::string> const& side_set_names)
{
  bool is_valid = true;

  dim_ = dim;
  std::string       bc_type_string("undefined");
  std::string       coordinate_string("undefined");
  std::stringstream ss(bc_string);
  ss >> bc_type_string;

  if (bc_type_string == "initial_velocity") {
    bc_type_ = INITIAL_VELOCITY;
  } else if (bc_type_string == "prescribed_velocity") {
    bc_type_ = PRESCRIBED_VELOCITY;
  } else if (bc_type_string == "prescribed_displacement") {
    bc_type_ = PRESCRIBED_DISPLACEMENT;
  } else if (bc_type_string == "prescribed_traction") {
    bc_type_ = PRESCRIBED_TRACTION;
  } else {
    throw std::logic_error(
        "Error processing boundary condition, unknown boundary condition "
        "type: " +
        bc_type_string);
  }

  bool const is_neumann_bc = bc_type_ == PRESCRIBED_TRACTION;

  if (is_neumann_bc == true) {
    ss >> side_set_name_;
  } else {
    ss >> node_set_name_;
  }
  ss >> coordinate_string;

  // figure out if magnitude is a double or an expression (check for quotes)
  int num_quotes = std::count(bc_string.begin(), bc_string.end(), '"');
  if (num_quotes == 2) {
    has_expression_    = true;
    size_t first_pos   = bc_string.find('"');
    size_t last_pos    = bc_string.rfind('"');
    expression_string_ = bc_string.substr(first_pos + 1, last_pos - first_pos - 1);
    expression_        = ExpressionParsing::BoundaryConditionFunctor(expression_string_);
  } else if (num_quotes == 0) {
    has_expression_ = false;
    ss >> magnitude_;  // DJL valgrind complains about this
  } else {
    throw std::logic_error("Error processing boundary condition, illegal number of quotes: " + bc_string);
  }

  std::transform(bc_type_string.begin(), bc_type_string.end(), bc_type_string.begin(), ::tolower);
  std::transform(coordinate_string.begin(), coordinate_string.end(), coordinate_string.begin(), ::tolower);

  // If id is -1, then either
  // 1) no nodes/sides on this processor belong to the node/side set, or
  // 2) the node/side set doesn't exist at all.
  if (is_neumann_bc == true) {
    side_set_id_ = find_name_id(side_set_names, side_set_name_);
    if (side_set_id_ == -1) is_valid = false;
  } else {
    node_set_id_ = find_name_id(node_set_names, node_set_name_);
    if (node_set_id_ == -1) is_valid = false;
  }

  if (coordinate_string == "x") {
    coordinate_ = 0;
  } else if (coordinate_string == "y") {
    coordinate_ = 1;
  } else if (coordinate_string == "z") {
    coordinate_ = 2;
  } else {
    throw std::logic_error("Error processing boundary condition, unknown coordinate: " + coordinate_string);
  }

  return is_valid;
}

}  // namespace nimble
