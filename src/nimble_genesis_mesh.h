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

#ifndef NIMBLE_GENESIS_MESH_H
#define NIMBLE_GENESIS_MESH_H

#ifdef NIMBLE_HAVE_DARMA
  #include "darma.h"
#else
  #include <map>
  #include <vector>
#endif

#include <string>


namespace nimble {

  class GenesisMesh {

  public:

    GenesisMesh() : file_name_("undefined"), dim_(-1) {}

    ~GenesisMesh() {}

#ifdef NIMBLE_HAVE_DARMA
    template<typename ArchiveType>
    void serialize(ArchiveType& ar) {
      ar | file_name_ | dim_ | node_global_id_ | node_x_ | node_y_ | node_z_ | elem_global_id_ | block_ids_ | block_names_ | block_elem_global_ids_ | block_num_nodes_per_elem_ | block_elem_connectivity_ | node_set_ids_ | node_set_names_ | node_sets_ | distribution_factors_;
    }
#endif

    bool IsValid() const {
      if (file_name_ == "none") {
        return false;
      }
      return true;
    }

    unsigned int GetNumNodes() const { return node_x_.size(); }

    const int * GetNodeGlobalIds() const { return &node_global_id_[0]; }

    unsigned int GetNumElements() const { return elem_global_id_.size(); }

    const int * GetElementGlobalIds() const { return &elem_global_id_[0]; }

    std::vector<int> const & GetElementGlobalIdsInBlock(int block_id) const { return block_elem_global_ids_.at(block_id); }

    unsigned int GetNumBlocks() const { return block_ids_.size(); }

    bool HasBlock(std::string const & block_name) const;

    std::vector<int> GetBlockIds() const { return block_ids_; }

    int GetNumElementsInBlock(int block_id) const;

    std::map<int, int> GetNumElementsInBlock() const;

    int GetNumNodesPerElement(int block_id) const {
      return block_num_nodes_per_elem_.at(block_id);
    }

    std::map<int, int> GetNumNodesPerElement() const { return block_num_nodes_per_elem_; }

    std::string GetElementType(int block_id) const;

    std::string GetBlockName(int block_id) const { return block_names_.at(block_id); }

    int GetBlockId(std::string const & block_name) const;

    void BlockNamesToOnProcessorBlockIds(std::vector<std::string> const & block_names,
                                         std::vector<int> & block_ids);

    int GetDim() const { return dim_; }

    const double * GetCoordinatesX() const { return &node_x_[0]; }

    const double * GetCoordinatesY() const { return &node_y_[0]; }

    const double * GetCoordinatesZ() const { return &node_z_[0]; }

    const int * GetConnectivity(int block_id) const { return &block_elem_connectivity_.at(block_id)[0]; }

    std::map<int, std::vector<int> > & GetConnectivity() { return block_elem_connectivity_; }

    int GetNumNodeSets() const { return static_cast<int>(node_set_ids_.size()); }

    std::vector<int> GetNodeSetIds() const { return node_set_ids_; }

    std::map<int, std::string> GetNodeSetNames() const { return node_set_names_; }

    std::map<int, std::vector<int> > GetNodeSets() const { return node_sets_; }

    std::map<int, std::vector<double> > GetDistributionFactors() const { return distribution_factors_; }

    void BoundingBox(double& x_min,
                     double& x_max,
                     double& y_min,
                     double& y_max,
                     double& z_min,
                     double& z_max) const ;

    std::vector<double> BoundingBoxCenter() const ;

    void AppendPeriodicPair(int local_master_node_id,
                            int local_slave_node_id,
                            const int * const global_node_ids,
                            std::map<int, int>& global_node_id_slave_to_master) const ;

    void CreatePeriodicRVELinearSystemMap(const int * const global_node_ids,
                                          std::vector<int>& linear_system_node_ids,
                                          std::map<int, std::vector<int>>& map_from_linear_system,
                                          int& corner_node_id) const ;

    void Print(bool verbose = false, int my_rank = 0) const ;

    void ReadFile(std::string file_name);

    void ReadTextFile(std::string file_name);

  protected:

    void ReportExodusError(int error_code, const char *method_name, const char *exodus_method_name) const;

    std::string file_name_;
    int dim_;

    std::vector<int> node_global_id_;
    std::vector<double> node_x_;
    std::vector<double> node_y_;
    std::vector<double> node_z_;
    std::vector<int> elem_global_id_;
    std::vector<int> block_ids_;
    std::map<int, std::string> block_names_;
    std::map<int, std::vector<int> > block_elem_global_ids_;
    std::map<int, int> block_num_nodes_per_elem_;
    std::map<int, std::vector<int> > block_elem_connectivity_;
    std::vector<int> node_set_ids_;
    std::map<int, std::string> node_set_names_;
    std::map<int, std::vector<int> > node_sets_;
    std::map<int, std::vector<double> > distribution_factors_;
  };

} // namespace nimble

#endif // NIMBLE_INPUT_EXODUS_H
