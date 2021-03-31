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

#ifndef NIMBLE_DATA_UTILS_H
#define NIMBLE_DATA_UTILS_H

#ifdef NIMBLE_HAVE_DARMA
#include "darma.h"
#else
#include <map>
#include <string>
#include <vector>
#endif

namespace nimble {

enum Relation
{
  UNDEFINED_RELATION = 0,
  NODE,
  ELEMENT,
  GLOBAL
};

enum Length
{
  LENGTH_0  = 0,
  LENGTH_1  = 1,
  LENGTH_2  = 2,
  LENGTH_3  = 3,
  LENGTH_4  = 4,
  LENGTH_5  = 5,
  LENTGH_6  = 6,
  LENGTH_7  = 7,
  LENGTH_8  = 8,
  LENGTH_9  = 9,
  LENGTH_10 = 10,
  LENGTH_11 = 11,
  LENGTH_12 = 12,
  LENGTH_13 = 13,
  LENGTH_14 = 14,
  LENGTH_15 = 15,
  LENGTH_16 = 16,
  LENGTH_17 = 17,
  LENGTH_18 = 18,
  LENGTH_19 = 19,
  LENGTH_20 = 20,
  SCALAR,
  VECTOR,
  SYMMETRIC_TENSOR,
  FULL_TENSOR,
  UNDEFINED_LENGTH
};

enum Step
{
  UNDEFINED_STEP = 0,
  STEP_N,
  STEP_NP1
};

int
LengthToInt(Length length, int dim);

std::string
AddIntegrationPointPrefix(std::string label, int ipt_number);

std::string
RemoveIntegrationPointPrefix(std::string label);

bool
HasIntegrationPointPrefix(std::string label);

int
LabelToIntegrationPointNumber(std::string label);

std::vector<std::string>
GetComponentLabels(std::string label, Length length, int dim);

std::string
GetComponentLabel(
    std::string label,
    Length      length,
    int         dim,
    int         component_index,
    int         ipt);

class Field
{
 public:
  Field()
      : id_(-1),
        label_("none"),
        relation_(UNDEFINED_RELATION),
        length_(UNDEFINED_LENGTH)
  {
  }

  virtual ~Field() {}

#ifdef NIMBLE_HAVE_DARMA
  template <typename ArchiveType>
  void
  serialize(ArchiveType& ar)
  {
    ar | id_ | label_ | relation_ | length_;
  }
#endif

  int         id_;
  std::string label_;
  Relation    relation_;
  Length      length_;
};

Length
LabelToLength(
    std::string                 label,
    std::map<int, Field> const& data_fields,
    int                         dim);

}  // namespace nimble

#endif
