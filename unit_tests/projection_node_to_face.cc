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

/*
 * projection_node_to_face.cc
 */

#include <iostream>
#include <string>

#include <gtest/gtest.h>

#include <nimble_contact_entity.h>
#include <nimble_contact_manager.h>


TEST(UnitTest, Projection) {

  using namespace nimble;

  double nCoord[3] = {0.123, 0.456, 0.789};
  ContactEntity node(ContactEntity::NODE, 10, nCoord, 0.0001, 0);

  double coord[9] = {0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0};
  int dummy_map[4] = {2, 2, 2, 2};

  ContactEntity face(ContactEntity::TRIANGLE, 0, coord, 1.0, 0, 1, dummy_map);

  ContactManager::PROJECTION_TYPE flag = ContactManager::PROJECTION_TYPE::UNKNOWN;
  ContactEntity::vertex projectedPoint;
  double gap = 0.0;
  double normal[3] = {0.0, 0.0, 0.0};

  double tol = 1.0e-08;
  ContactManager::SimpleClosestPointProjectionSingle(node, face,
     &flag, &projectedPoint, gap, normal, tol);

  EXPECT_EQ(flag, ContactManager::PROJECTION_TYPE::FACE);

  EXPECT_DOUBLE_EQ(projectedPoint.coords_[0], nCoord[0]);
  EXPECT_DOUBLE_EQ(projectedPoint.coords_[1], nCoord[1]);
  EXPECT_DOUBLE_EQ(projectedPoint.coords_[2], 0.0);

  EXPECT_DOUBLE_EQ(gap, nCoord[2]);

  //----

  nCoord[0] = 1.0; nCoord[1] = 0.0; nCoord[2] = -1.23;
  node.SetCoordinates(nCoord);
  ContactManager::SimpleClosestPointProjectionSingle(node, face,
     &flag, &projectedPoint, gap, normal, tol);

  EXPECT_EQ(flag, ContactManager::PROJECTION_TYPE::FACE);

  EXPECT_DOUBLE_EQ(projectedPoint.coords_[0], nCoord[0]);
  EXPECT_DOUBLE_EQ(projectedPoint.coords_[1], nCoord[1]);
  EXPECT_DOUBLE_EQ(projectedPoint.coords_[2], 0.0);

  EXPECT_DOUBLE_EQ(gap, nCoord[2]);

  //----

  nCoord[0] = 1.0; nCoord[1] = 0.5 * tol; nCoord[2] = 0.777;
  node.SetCoordinates(nCoord);
  ContactManager::SimpleClosestPointProjectionSingle(node, face,
     &flag, &projectedPoint, gap, normal, tol);

  EXPECT_EQ(flag, ContactManager::PROJECTION_TYPE::FACE);

  EXPECT_DOUBLE_EQ(projectedPoint.coords_[0], nCoord[0]);
  EXPECT_DOUBLE_EQ(projectedPoint.coords_[1], nCoord[1]);
  EXPECT_DOUBLE_EQ(projectedPoint.coords_[2], 0.0);

  EXPECT_DOUBLE_EQ(gap, nCoord[2]);

  //----

  nCoord[0] = 1.0; nCoord[1] = 1.11; nCoord[2] = -1.23;
  node.SetCoordinates(nCoord);
  ContactManager::SimpleClosestPointProjectionSingle(node, face,
     &flag, &projectedPoint, gap, normal, tol);

  EXPECT_EQ(flag, ContactManager::PROJECTION_TYPE::UNKNOWN);

}

