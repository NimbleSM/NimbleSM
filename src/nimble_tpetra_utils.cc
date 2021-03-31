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

#include "nimble_tpetra_utils.h"

#include <set>

#include "nimble_genesis_mesh.h"
#include "nimble_model_data.h"

namespace nimble {

void
TpetraContainer::Initialize(
    GenesisMesh const&      mesh,
    std::vector<int> const& global_node_ids,
    comm_type               comm)
{
  Initialize(mesh.GetDim(), mesh.GetNumNodes(), comm, global_node_ids);
}

void
TpetraContainer::Initialize(
    int                     d,
    unsigned int            n,
    comm_type               comm,
    std::vector<int> const& global_node_ids)
{
#ifdef NIMBLE_HAVE_TRILINOS

  comm_            = comm;
  global_node_ids_ = global_node_ids;

  unsigned int num_nodes = n;
  dim_                   = d;

  const Tpetra::global_size_t num_global_elements =
      Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
  const global_ordinal_type index_base = 0;

  SimpleTieBreak<local_ordinal_type, global_ordinal_type> tie_break;

  // tpetra vector for scalar data (e.g., lumped mass)
  local_ordinal_type               num_1d_entries = num_nodes;
  std::vector<global_ordinal_type> my_1d_entries(num_1d_entries);
  for (int i_node = 0; i_node < num_nodes; i_node++) {
    my_1d_entries[i_node] = global_node_ids_[i_node];
  }
  Teuchos::RCP<const map_type> map_1d = Teuchos::rcp(new map_type(
      num_global_elements,
      my_1d_entries.data(),
      num_1d_entries,
      index_base,
      comm_));
  vec_1d_                             = Teuchos::rcp(new vector_type(map_1d));

  Teuchos::RCP<const map_type> map_1d_one_to_one =
      Tpetra::createOneToOne(map_1d, tie_break);
  vec_1d_one_to_one_ = Teuchos::rcp(new vector_type(map_1d_one_to_one));

  export_into_vec_1d_one_to_one_ =
      Teuchos::rcp(new export_type(map_1d, map_1d_one_to_one));
  import_from_vec_1d_one_to_one_ =
      Teuchos::rcp(new import_type(map_1d_one_to_one, map_1d));

  // tpetra vector for vector data (e.g., internal force)
  local_ordinal_type               num_3d_entries = num_nodes * dim_;
  std::vector<global_ordinal_type> my_3d_entries(num_3d_entries);
  for (int i_node = 0; i_node < num_nodes; i_node++) {
    for (int dof = 0; dof < dim_; dof++) {
      my_3d_entries[i_node * dim_ + dof] =
          global_node_ids_[i_node] * dim_ + dof;
    }
  }
  Teuchos::RCP<const map_type> map_3d = Teuchos::rcp(new map_type(
      num_global_elements,
      my_3d_entries.data(),
      num_3d_entries,
      index_base,
      comm_));
  vec_3d_                             = Teuchos::rcp(new vector_type(map_3d));

  Teuchos::RCP<const map_type> map_3d_one_to_one =
      Tpetra::createOneToOne(map_3d, tie_break);
  vec_3d_one_to_one_ = Teuchos::rcp(new vector_type(map_3d_one_to_one));

  export_into_vec_3d_one_to_one_ =
      Teuchos::rcp(new export_type(map_3d, map_3d_one_to_one));
  import_from_vec_3d_one_to_one_ =
      Teuchos::rcp(new import_type(map_3d_one_to_one, map_3d));

#endif
}

void
TpetraContainer::AllocateTangentStiffnessMatrix(GenesisMesh const& mesh)
{
#ifdef NIMBLE_HAVE_TRILINOS

  // record the nonzero entries in each row
  // the map key is the global id of the row
  // the entries in the set are the global ids of the entries in a given row
  using global_ordinal_set = std::set<global_ordinal_type>;
  std::map<global_ordinal_type, global_ordinal_set> matrix_nonzeros;

  int              dim       = mesh.GetDim();
  std::vector<int> block_ids = mesh.GetBlockIds();

  for (auto const& block_id : block_ids) {
    int              num_elem_in_block  = mesh.GetNumElementsInBlock(block_id);
    int              num_nodes_per_elem = mesh.GetNumNodesPerElement(block_id);
    const int* const conn               = mesh.GetConnectivity(block_id);
    int              conn_offset        = 0;
    for (int i_elem = 0; i_elem < num_elem_in_block; i_elem++) {
      for (int i_node = 0; i_node < num_nodes_per_elem; i_node++) {
        global_ordinal_type global_id_i =
            global_node_ids_[conn[conn_offset + i_node]];
        for (int j_node = 0; j_node < num_nodes_per_elem; j_node++) {
          global_ordinal_type global_id_j =
              global_node_ids_[conn[conn_offset + j_node]];
          for (int i_dof = 0; i_dof < dim; i_dof++) {
            for (int j_dof = 0; j_dof < dim; j_dof++) {
              global_ordinal_type matrix_entry_i = global_id_i * dim + i_dof;
              global_ordinal_type matrix_entry_j = global_id_j * dim + j_dof;
              auto&& row_nonzeros = matrix_nonzeros[matrix_entry_i];
              row_nonzeros.insert(matrix_entry_j);
            }
          }
        }
      }
      conn_offset += num_nodes_per_elem;
    }
  }

  // record the rows that have entries on this processor and count nonzeros
  std::vector<global_ordinal_type> my_entries(matrix_nonzeros.size());
  {
    int row_index = 0;
    for (auto const& entry : matrix_nonzeros) {
      my_entries[row_index] = entry.first;
      ++row_index;
    }
  }
  std::sort(my_entries.begin(), my_entries.end());
  std::vector<global_ordinal_set::size_type> num_matrix_nonzeros_per_row(
      matrix_nonzeros.size());
  {
    int row_index = 0;
    for (auto&& entry : my_entries) {
      num_matrix_nonzeros_per_row[row_index++] =
          matrix_nonzeros.at(entry).size();
    }
  }

  // create a row map
  const Tpetra::global_size_t num_global_elements =
      Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
  const global_ordinal_type    index_base = 0;
  Teuchos::RCP<const map_type> row_map    = Teuchos::rcp(new map_type(
      num_global_elements,
      my_entries.data(),
      static_cast<local_ordinal_type>(my_entries.size()),
      index_base,
      comm_));

  // create a crs graph based on the row map
  Teuchos::RCP<graph_type> crs_graph = Teuchos::rcp(new graph_type(
      row_map, Teuchos::arrayViewFromVector(num_matrix_nonzeros_per_row)));

  // identify the nonzeros in the crs graph
  for (auto const& entry : matrix_nonzeros) {
    auto                             global_row          = entry.first;
    auto&&                           nonzeros_in_row_set = entry.second;
    std::vector<global_ordinal_type> nonzeros_in_row(
        nonzeros_in_row_set.begin(), nonzeros_in_row_set.end());
    crs_graph->insertGlobalIndices(
        global_row, Teuchos::arrayViewFromVector(nonzeros_in_row));
  }
  crs_graph->fillComplete();

  // allocate the tangent stiffness matrix
  crs_matrix_           = Teuchos::rcp(new matrix_type(crs_graph));
  crs_matrix_container_ = TpetraMatrixContainer<
      scalar_type,
      local_ordinal_type,
      global_ordinal_type>(crs_matrix_);

#endif
}

int
TpetraContainer::TangentStiffnessMatrixNumNonzeros()
{
#ifdef NIMBLE_HAVE_TRILINOS
  return crs_matrix_->getGlobalNumEntries();
#else
  return 0;
#endif
}

void
TpetraContainer::TangentStiffnessMatrixSetScalar(double value)
{
#ifdef NIMBLE_HAVE_TRILINOS
  crs_matrix_->setAllToScalar(value);
#endif
}

void
TpetraContainer::TangentStiffnessMatrixReplaceValue(
    int    row,
    int    col,
    double value)
{
#ifdef NIMBLE_HAVE_TRILINOS

  local_ordinal_type local_row = static_cast<local_ordinal_type>(row);
  Teuchos::ArrayRCP<local_ordinal_type> local_cols;
  Teuchos::ArrayRCP<scalar_type>        values;

  local_ordinal_type success =
      crs_matrix_->sumIntoLocalValues(local_row, local_cols(), values());

  if (success == Teuchos::OrdinalTraits<local_ordinal_type>::invalid()) {
    throw std::logic_error(
        "\nError in TangentStiffMatrixReplaceValue(), replaceLocalValues() "
        "failed.\n");
  } else if (success != 1) {
    throw std::logic_error(
        "\nError in TangentStiffMatrixReplaceValue(), invalid index.\n");
  }
  std::cout << "DJL DEBUGGING " << success << std::endl;

#endif
}

void
TpetraContainer::VectorReduction(
    ModelData&  model_data,
    std::string quantity_label)
{
#ifdef NIMBLE_HAVE_TRILINOS
  int     field_id = model_data.GetFieldId(quantity_label);
  Field   field    = model_data.GetField(field_id);
  Length  length   = field.length_;
  double* data     = model_data.GetNodeData(field_id);
  //
  int data_dimension = LengthToInt(length, dim_);
  VectorReduction(data_dimension, data);
#endif
}

void
TpetraContainer::VectorReduction(int data_dimension, double* data)
{
#ifdef NIMBLE_HAVE_TRILINOS
  Length length = SCALAR;
  if (data_dimension == dim_) length = VECTOR;

  if (length == SCALAR) {
    // load the data into a tpetra vector
    Teuchos::ArrayRCP<scalar_type> vec_1d_data = vec_1d_->getDataNonConst();
    for (unsigned int i_node = 0; i_node < global_node_ids_.size(); i_node++) {
      vec_1d_data[i_node] = data[i_node];
    }

    // use tpetra import/export to allreduce over the vector
    vec_1d_one_to_one_->doExport(
        *vec_1d_, *export_into_vec_1d_one_to_one_, Tpetra::CombineMode::ADD);
    vec_1d_->doImport(
        *vec_1d_one_to_one_,
        *import_from_vec_1d_one_to_one_,
        Tpetra::CombineMode::INSERT);

    // copy the reduced values back into the original container
    for (unsigned int i_node = 0; i_node < global_node_ids_.size(); i_node++) {
      data[i_node] = vec_1d_data[i_node];
    }
  } else if (length == VECTOR) {
    // load the data into a tpetra vector
    Teuchos::ArrayRCP<scalar_type> vec_3d_data = vec_3d_->getDataNonConst();
    int                            num_dof     = LengthToInt(length, dim_);
    for (unsigned int i_node = 0; i_node < global_node_ids_.size(); i_node++) {
      for (int dof = 0; dof < num_dof; dof++) {
        vec_3d_data[i_node * dim_ + dof] = data[i_node * dim_ + dof];
      }
    }

    // use tpetra import/export to allreduce over the vector
    vec_3d_one_to_one_->doExport(
        *vec_3d_, *export_into_vec_3d_one_to_one_, Tpetra::CombineMode::ADD);
    vec_3d_->doImport(
        *vec_3d_one_to_one_,
        *import_from_vec_3d_one_to_one_,
        Tpetra::CombineMode::INSERT);

    // copy the reduced values back into the original container
    for (unsigned int i_node = 0; i_node < global_node_ids_.size(); i_node++) {
      for (int dof = 0; dof < num_dof; dof++) {
        data[i_node * dim_ + dof] = vec_3d_data[i_node * dim_ + dof];
      }
    }
  }

#endif
}

}  // namespace nimble
