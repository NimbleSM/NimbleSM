
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

#include "nimble_contact_entity.h"

namespace nimble {

void
ContactEntity::ExportGeometryInto(ContactEntity& xerox) const
{
  xerox.entity_type_ = entity_type_;

  xerox.num_nodes_ = num_nodes_;
  xerox.coord_1_x_ = coord_1_x_;
  xerox.coord_1_y_ = coord_1_y_;
  xerox.coord_1_z_ = coord_1_z_;
  xerox.coord_2_x_ = coord_2_x_;
  xerox.coord_2_y_ = coord_2_y_;
  xerox.coord_2_z_ = coord_2_z_;
  xerox.coord_3_x_ = coord_3_x_;
  xerox.coord_3_y_ = coord_3_y_;
  xerox.coord_3_z_ = coord_3_z_;

  xerox.char_len_ = char_len_;

  xerox.bounding_box_x_min_ = bounding_box_x_min_;
  xerox.bounding_box_x_max_ = bounding_box_x_max_;
  xerox.bounding_box_y_min_ = bounding_box_y_min_;
  xerox.bounding_box_y_max_ = bounding_box_y_max_;
  xerox.bounding_box_z_min_ = bounding_box_z_min_;
  xerox.bounding_box_z_max_ = bounding_box_z_max_;
}

void
ContactEntity::SetBoundingBox()
{
  centroid_[0]        = coord_1_x_;
  centroid_[1]        = coord_1_y_;
  centroid_[2]        = coord_1_z_;
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
    if (coord_2_x_ < bounding_box_x_min_) bounding_box_x_min_ = coord_2_x_;
    if (coord_2_x_ > bounding_box_x_max_) bounding_box_x_max_ = coord_2_x_;
    if (coord_2_y_ < bounding_box_y_min_) bounding_box_y_min_ = coord_2_y_;
    if (coord_2_y_ > bounding_box_y_max_) bounding_box_y_max_ = coord_2_y_;
    if (coord_2_z_ < bounding_box_z_min_) bounding_box_z_min_ = coord_2_z_;
    if (coord_2_z_ > bounding_box_z_max_) bounding_box_z_max_ = coord_2_z_;
    centroid_[0] += coord_3_x_;
    centroid_[1] += coord_3_y_;
    centroid_[2] += coord_3_z_;
    if (coord_3_x_ < bounding_box_x_min_) bounding_box_x_min_ = coord_3_x_;
    if (coord_3_x_ > bounding_box_x_max_) bounding_box_x_max_ = coord_3_x_;
    if (coord_3_y_ < bounding_box_y_min_) bounding_box_y_min_ = coord_3_y_;
    if (coord_3_y_ > bounding_box_y_max_) bounding_box_y_max_ = coord_3_y_;
    if (coord_3_z_ < bounding_box_z_min_) bounding_box_z_min_ = coord_3_z_;
    if (coord_3_z_ > bounding_box_z_max_) bounding_box_z_max_ = coord_3_z_;
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

#ifdef NIMBLE_HAVE_BVH
// Free functions for accessing entity info for bvh
bvh::bphase_kdop
get_entity_kdop(const ContactEntity& _entity)
{
  return _entity.Kdop();
}

std::size_t
get_entity_global_id(const ContactEntity& _entity)
{
  return static_cast<std::size_t>(_entity.contact_entity_global_id());
}

bvh::m::vec3d
get_entity_centroid(const ContactEntity& _entity)
{
  const auto centroid = _entity.centroid();

  return bvh::m::vec3d{centroid[0], centroid[1], centroid[2]};
}
#endif

}  // namespace nimble
