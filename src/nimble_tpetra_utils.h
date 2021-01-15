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

#ifndef NIMBLE_TPETRA_UTILS_H
#define NIMBLE_TPETRA_UTILS_H

#include "nimble_genesis_mesh.h"
#include "nimble_data_manager.h"

#ifdef NIMBLE_HAVE_TRILINOS
  #include <Tpetra_Core.hpp>
  #include <Teuchos_GlobalMPISession.hpp>
  #include <Teuchos_oblackholestream.hpp>
  #include <Tpetra_Vector.hpp>
  #include <Tpetra_CrsMatrix.hpp>
  #include <Tpetra_TieBreak.hpp>
#endif

// TODO
// Move the matrix contain class to nimble_linear_solver.h
// Modify accessor routines for all matrix container classes to be as compatible as possible with CrsMatrix
// Implement row-wise inserts, to be as efficient as possible for the epetra case (probably all cases)

namespace nimble {

#ifdef NIMBLE_HAVE_TRILINOS
  typedef Teuchos::RCP<const Teuchos::Comm<int>> comm_type;
#else
  typedef int comm_type;
#endif

  template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal>
  class TpetraMatrixContainer {

   public:

    TpetraMatrixContainer() {}

    ~TpetraMatrixContainer() {}

#ifdef NIMBLE_HAVE_TRILINOS

    TpetraMatrixContainer(Teuchos::RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal>> crs_matrix)
      : crs_matrix_(crs_matrix) {}

    bool SumIntoValue(LocalOrdinal local_row_index, GlobalOrdinal local_col_index, Scalar value) {

      LocalOrdinal local_row = static_cast<LocalOrdinal>(local_row_index);
      Teuchos::ArrayRCP<LocalOrdinal> local_cols;
      Teuchos::ArrayRCP<Scalar> values;

      LocalOrdinal success = crs_matrix_->sumIntoLocalValues(local_row, local_cols(), values());

      if (success == Teuchos::OrdinalTraits<LocalOrdinal>::invalid()) {
        throw std::logic_error("\nError in TpetraMatrixContainer::sumIntoValue(), replaceLocalValues() failed.\n");
      }
      else if (success != 1) {
        throw std::logic_error("\nError in TpetraMatrixContainer::sumIntoValue(), invalid index.\n");
      }
      std::cout << "DJL DEBUGGING TpetraMatrixContainer::sumIntoValue() success! " << success << std::endl;
    }

  private:

    Teuchos::RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal>> crs_matrix_;

#endif
  };

  class TpetraContainer {

  public:

#ifdef NIMBLE_HAVE_TRILINOS
    typedef Tpetra::Map<> map_type;
    typedef Tpetra::Vector<> vector_type;
    typedef Tpetra::Vector<>::scalar_type scalar_type;
    typedef Tpetra::Vector<>::local_ordinal_type local_ordinal_type;
    typedef Tpetra::Vector<>::global_ordinal_type global_ordinal_type;
    typedef Tpetra::Import<local_ordinal_type, global_ordinal_type> import_type;
    typedef Tpetra::Export<local_ordinal_type, global_ordinal_type> export_type;
    typedef Tpetra::CrsGraph<local_ordinal_type, global_ordinal_type> graph_type;
    typedef Tpetra::CrsMatrix<scalar_type, local_ordinal_type, global_ordinal_type> matrix_type;

    template <typename LocalOrdinal, typename GlobalOrdinal>
    class SimpleTieBreak : public Tpetra::Details::TieBreak<LocalOrdinal, GlobalOrdinal> {
      std::size_t
      selectedIndex (GlobalOrdinal GID,
                     const std::vector<std::pair<int, LocalOrdinal> >& pid_and_lid) const {
        std::size_t selected_index = 0;
        int selected_pid = pid_and_lid[0].first;
        for (unsigned int i=0 ; i<pid_and_lid.size() ; i++) {
          int pid = pid_and_lid[i].first;
          if (pid < selected_pid) {
            selected_index = i;
            selected_pid = pid;
          }
        }
        return selected_index;
      }
    };
#endif

    TpetraContainer() {}

    void Initialize(GenesisMesh const & mesh,
                    std::vector<int> const & global_node_ids,
                    comm_type comm);

  void Initialize(int d, unsigned int n, comm_type comm,
                  std::vector<int> const & global_node_ids);

  void AllocateTangentStiffnessMatrix(GenesisMesh const & mesh);

    int TangentStiffnessMatrixNumNonzeros();

    void TangentStiffnessMatrixSetScalar(double value);

    void TangentStiffnessMatrixReplaceValue(int row, int col, double value);

#ifdef NIMBLE_HAVE_TRILINOS
    TpetraMatrixContainer<scalar_type, local_ordinal_type, global_ordinal_type>& GetTpetraMatrixContainer() {
      return crs_matrix_container_;
    }
#endif

    void VectorReduction(ModelData & model_data,
                         std::string quantity_label);

  void VectorReduction(int data_dimension, double* data);

  int dim_;
    std::vector<int> global_node_ids_;

#ifdef NIMBLE_HAVE_TRILINOS
    comm_type comm_;
    Teuchos::RCP<vector_type> vec_1d_;
    Teuchos::RCP<vector_type> vec_1d_one_to_one_;
    Teuchos::RCP<vector_type> vec_3d_;
    Teuchos::RCP<vector_type> vec_3d_one_to_one_;
    Teuchos::RCP<export_type> export_into_vec_1d_one_to_one_;
    Teuchos::RCP<import_type> import_from_vec_1d_one_to_one_;
    Teuchos::RCP<export_type> export_into_vec_3d_one_to_one_;
    Teuchos::RCP<import_type> import_from_vec_3d_one_to_one_;
    Teuchos::RCP<matrix_type> crs_matrix_;
    TpetraMatrixContainer<scalar_type, local_ordinal_type, global_ordinal_type> crs_matrix_container_;
#endif
  };

} // namespace nimble

#endif
