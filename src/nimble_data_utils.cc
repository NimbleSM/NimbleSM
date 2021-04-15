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

#include "nimble_data_utils.h"

#include <ctype.h>

#include <algorithm>
#include <sstream>
#include <stdexcept>

namespace nimble {

int
LengthToInt(Length length, int dim)
{
  int value;
  if (dim == 2) {
    switch (length) {
      case SCALAR: value = 1; break;
      case VECTOR: value = 2; break;
      case SYMMETRIC_TENSOR: value = 3; break;
      case FULL_TENSOR: value = 4; break;
      case UNDEFINED_LENGTH: throw std::logic_error("\nError in LengthToInt(), unrecognized Length.\n");
      default: value = static_cast<int>(length); break;
    }
  } else if (dim == 3) {
    switch (length) {
      case SCALAR: value = 1; break;
      case VECTOR: value = 3; break;
      case SYMMETRIC_TENSOR: value = 6; break;
      case FULL_TENSOR: value = 9; break;
      case UNDEFINED_LENGTH: throw std::logic_error("\nError in LengthToInt(), unrecognized Length.\n");
      default: value = static_cast<int>(length); break;
    }
  } else {
    throw std::logic_error("\nError in LengthToInt(), unrecognized dimension.\n");
  }
  return value;
}

std::string
AddIntegrationPointPrefix(std::string label, int ipt_number)
{
  std::string       modified_label;
  std::stringstream ss;
  ss << ipt_number;
  if (ss.str().size() == 1) {
    modified_label = "ipt0" + ss.str() + "_" + label;
  } else {
    modified_label = "ipt" + ss.str() + "_" + label;
  }
  return modified_label;
}

std::string
RemoveIntegrationPointPrefix(std::string label)
{
  std::string stripped_label(label);

  size_t pos = label.find("ipt");
  if (pos == 0) {
    // remove "ipt"
    std::string temp = label.substr(3);
    // strip off integers, they are the integration point number
    // the initial value of 1 is for the trailing underscore
    int  num_digits   = 1;
    bool first_digits = true;
    for (char& c : temp) {
      if (first_digits && isdigit(c)) {
        num_digits++;
      } else {
        first_digits = false;
      }
    }
    stripped_label = temp.substr(num_digits);
  }
  return stripped_label;
}

bool
HasIntegrationPointPrefix(std::string label)
{
  std::string label_without_prefix = RemoveIntegrationPointPrefix(label);
  if (label == label_without_prefix) { return false; }
  return true;
}

int
LabelToIntegrationPointNumber(std::string label)
{
  int int_pt_number;

  size_t pos = label.find("ipt");
  if (pos != 0) {
    throw std::logic_error(
        "\nError:  LabelToIntegratinPointNumber() called with a label that "
        "does not beging with ipt.\n");
  }

  // remove "ipt"
  std::string temp = label.substr(3);

  std::stringstream ss;
  bool              first_digits = true;
  for (char& c : temp) {
    if (first_digits && isdigit(c)) {
      ss << c;
    } else {
      first_digits = false;
    }
  }
  ss >> int_pt_number;

  return int_pt_number;
}

std::vector<std::string>
GetComponentLabels(std::string label, Length length, int dim)
{
  std::vector<std::string> labels;
  if (length == SCALAR) {
    labels.push_back(label);
  } else if (length == VECTOR) {
    labels.push_back(label + "_x");
    labels.push_back(label + "_y");
    if (dim == 3) { labels.push_back(label + "_z"); }
  } else if (length == SYMMETRIC_TENSOR) {
    if (dim == 2) {
      labels.push_back(label + "_xx");
      labels.push_back(label + "_yy");
      labels.push_back(label + "_xy");
    } else if (dim == 3) {
      labels.push_back(label + "_xx");
      labels.push_back(label + "_yy");
      labels.push_back(label + "_zz");
      labels.push_back(label + "_xy");
      labels.push_back(label + "_yz");
      labels.push_back(label + "_zx");
    }
  } else if (length == FULL_TENSOR) {
    if (dim == 2) {
      labels.push_back(label + "_xx");
      labels.push_back(label + "_yy");
      labels.push_back(label + "_xy");
      labels.push_back(label + "_xy");
    } else if (dim == 3) {
      labels.push_back(label + "_xx");
      labels.push_back(label + "_yy");
      labels.push_back(label + "_zz");
      labels.push_back(label + "_xy");
      labels.push_back(label + "_yz");
      labels.push_back(label + "_zx");
      labels.push_back(label + "_yx");
      labels.push_back(label + "_zy");
      labels.push_back(label + "_xz");
    }
  } else {
    int num_labels = static_cast<int>(length);
    for (int i = 0; i < num_labels; i++) {
      std::stringstream ss;
      ss << "_" << i + 1;
      labels.push_back(label + ss.str());
    }
  }
  return labels;
}

std::string
GetComponentLabel(std::string label, Length length, int dim, int component_index, int ipt)
{
  if (length == VECTOR) {
    switch (component_index) {
      case 0: label += "_x"; break;
      case 1: label += "_y"; break;
      case 2: label += "_z"; break;
      default: throw std::logic_error("\nError in GetComponentLabel().\n");
    }
  } else if (length == SYMMETRIC_TENSOR) {
    if (dim == 2) {
      switch (component_index) {
        case 0: label += "_xx"; break;
        case 1: label += "_yy"; break;
        case 2: label += "_xy"; break;
        default: throw std::logic_error("\nError in GetComponentLabel().\n");
      }
    } else if (dim == 3) {
      switch (component_index) {
        case 0: label += "_xx"; break;
        case 1: label += "_yy"; break;
        case 2: label += "_zz"; break;
        case 3: label += "_xy"; break;
        case 4: label += "_yz"; break;
        case 5: label += "_zx"; break;
        default: throw std::logic_error("\nError in GetComponentLabel().\n");
      }
    }
  } else if (length == FULL_TENSOR) {
    if (dim == 2) {
      switch (component_index) {
        case 0: label += "_xx"; break;
        case 1: label += "_yy"; break;
        case 2: label += "_xy"; break;
        case 3: label += "_yx"; break;
        default: throw std::logic_error("\nError in GetComponentLabel().\n");
      }
    } else if (dim == 3) {
      switch (component_index) {
        case 0: label += "_xx"; break;
        case 1: label += "_yy"; break;
        case 2: label += "_zz"; break;
        case 3: label += "_xy"; break;
        case 4: label += "_yz"; break;
        case 5: label += "_zx"; break;
        case 6: label += "_yx"; break;
        case 7: label += "_zy"; break;
        case 8: label += "_xz"; break;
        default: throw std::logic_error("\nError in GetComponentLabel().\n");
      }
    }
  }

  label = AddIntegrationPointPrefix(label, ipt);

  return label;
}

Length
LabelToLength(std::string label, std::map<int, Field> const& data_fields, int dim)
{
  label = RemoveIntegrationPointPrefix(label);
  for (auto const& it : data_fields) {
    Field const& field = it.second;
    if (label == RemoveIntegrationPointPrefix(field.label_)) { return field.length_; }
  }
  throw std::logic_error("\nError in LabelToField(), failed to find label " + label + ".\n");
}

}  // namespace nimble
