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

#include "nimble_mesh_utils.h"

namespace nimble {

  void DetermineTangentMatrixNonzeroStructure(GenesisMesh const & mesh,
                                              std::vector<int> const & linear_system_node_ids,
                                              std::vector<int> & i_index,
                                              std::vector<int> & j_index) {

    std::vector<int> block_ids = mesh.GetBlockIds();
    int dim = mesh.GetDim();
    int num_rows = linear_system_node_ids.size() * dim;

    // column indexes for the nonzeros in each row
    std::vector< std::set<int> > nonzeros(num_rows);

    for (auto const & block_id : block_ids) {
      int num_elem = mesh.GetNumElementsInBlock(block_id);
      int num_nodes_per_elem = mesh.GetNumNodesPerElement(block_id);
      std::vector<int> linear_system_dof(num_nodes_per_elem * dim);
      int const * elem_conn = mesh.GetConnectivity(block_id);
      for (int i_elem=0 ; i_elem < num_elem ; i_elem++) {
        for (int i_node=0 ; i_node<num_nodes_per_elem ; i_node++) {
          for (int i_dim=0 ; i_dim<dim ; i_dim++) {
            linear_system_dof[i_node * dim + i_dim] = linear_system_node_ids[ elem_conn[num_nodes_per_elem*i_elem + i_node] ] * dim + i_dim;
          }
        }
        int row, col;
        for (int i=0 ; i<num_nodes_per_elem * dim ; i++) {
          for (int j=0 ; j<num_nodes_per_elem * dim ; j++) {
            row = linear_system_dof[i];
            col = linear_system_dof[j];
            nonzeros[row].insert(col);
            nonzeros[col].insert(row);
          }
        }
      }
    }

    // i_index and j_index are arrays containing the row and column indices, respectively, for each nonzero
    int num_entries(0);
    for (int i_row=0 ; i_row<num_rows ; i_row++) {
      num_entries += nonzeros[i_row].size();
    }
    i_index.resize(num_entries);
    j_index.resize(num_entries);
    int index(0);
    for (int i_row=0 ; i_row<num_rows ; i_row++) {
      for (auto const & entry : nonzeros[i_row]) {
        i_index[index] = i_row;
        j_index[index++] = entry;
      }
    }
  }

} // namespace nimble
