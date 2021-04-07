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

#ifndef NIMBLE_PARSER_H
#define NIMBLE_PARSER_H

#include <sstream>
#include <stdexcept>

#ifdef NIMBLE_HAVE_DARMA
#include "darma.h"
#else
#include <map>
#include <string>
#include <vector>
#endif

#ifdef NIMBLE_HAVE_TRILINOS
#include <Tpetra_Core.hpp>
#endif

namespace nimble {

constexpr int MAX_C_STR_SIZE = 256;

std::string
IOFileName(
    std::string const& serial_name,
    std::string const& extension,
    std::string const& label     = std::string(),
    int                my_rank   = 0,
    int                num_ranks = 0);

void
IOFileNameThreadSafe(
    char const* const serial_file_name,
    char const* const extension,
    char const* const label,
    int               my_mpi_rank,
    int               mpi_num_ranks,
    int               my_thread_rank,
    int               num_threads,
    char* const       file_name);

struct BlockProperties
{
  BlockProperties() : block_name_("none"), block_id_(-1) {}
  BlockProperties(std::string props);

#ifdef NIMBLE_HAVE_DARMA
  template <typename ArchiveType>
  void
  serialize(ArchiveType& ar)
  {
    ar | block_name_ | block_id_ | material_key_;
  }
#endif

  std::string block_name_;
  int         block_id_;
  std::string material_key_;
};

class Parser
{
 public:
  Parser();

  virtual ~Parser() = default;

#ifdef NIMBLE_HAVE_DARMA
  template <typename ArchiveType>
  void
  serialize(ArchiveType& ar)
  {
    ar | file_name_ | genesis_file_name_ | rve_genesis_file_name_;
    ar | exodus_file_name_ | use_two_level_mesh_decomposition_;
    ar | write_timing_data_file_ | time_integration_scheme_;
    ar | nonlinear_solver_relative_tolerance_ | nonlinear_solver_max_iterations_;
    ar | initial_time_ | final_time_ | num_load_steps_ | output_frequency_ | reduction_version_;
    ar | contact_string_ | visualize_contact_entities_ | visualize_contact_bounding_boxes_;
    ar | contact_visualization_file_name_ | microscale_output_element_ids_ | material_strings_;
    ar | macroscale_blocks_ | microscale_blocks_ | microscale_boundary_condition_strategy_;
    ar | boundary_condition_strings_ | output_field_string_;
    ar | file_name_;
    ar | env_set_ | use_kokkos_ | use_tpetra_ | use_vt_;
    ar | my_rank_ | num_ranks_;
  }
#endif

  void
  Initialize();

  std::string
  GenesisFileName() const
  {
    return genesis_file_name_;
  }

  std::string
  RVEGenesisFileName() const
  {
    return rve_genesis_file_name_;
  }

  std::string
  ExodusFileName() const
  {
    return exodus_file_name_;
  }

  bool
  UseTwoLevelMeshDecomposition() const
  {
    return use_two_level_mesh_decomposition_;
  }

  bool
  WriteTimingDataFile() const
  {
    return write_timing_data_file_;
  }

  std::string
  TimeIntegrationScheme() const
  {
    if (time_integration_scheme_ != "explicit" && time_integration_scheme_ != "quasistatic") {
      std::string msg =
          "\n**** Error in Parser::TimeIntegrationScheme(), invalid "
          "integration scheme " +
          time_integration_scheme_ + ".\n";
      throw std::logic_error(msg);
    }
    return time_integration_scheme_;
  }

  double
  NonlinearSolverRelativeTolerance() const
  {
    return nonlinear_solver_relative_tolerance_;
  }

  int
  NonlinearSolverMaxIterations() const
  {
    return nonlinear_solver_max_iterations_;
  }

  double
  InitialTime() const
  {
    return initial_time_;
  }

  double
  FinalTime() const
  {
    return final_time_;
  }

  int
  NumLoadSteps() const
  {
    return num_load_steps_;
  }

  int
  OutputFrequency() const
  {
    return output_frequency_;
  }

  bool
  HasContact() const
  {
    return (!contact_string_.empty());
  }

  bool
  ContactVisualization() const
  {
    bool visualize_contact = visualize_contact_entities_ || visualize_contact_bounding_boxes_;
    return visualize_contact;
  }

  std::string
  ContactVisualizationFileName() const
  {
    return contact_visualization_file_name_;
  }

  std::string
  ContactString() const
  {
    return contact_string_;
  }

  std::vector<int>
  MicroscaleOutputElementIds() const
  {
    return microscale_output_element_ids_;
  }

  std::string
  GetMacroscaleMaterialParameters(int block_id) const
  {
    if (macroscale_blocks_.find(block_id) == macroscale_blocks_.end()) {
      std::string none_str("none");
      return none_str;
    }
    BlockProperties const& block          = macroscale_blocks_.at(block_id);
    std::string const&     material_key   = block.material_key_;
    std::string const&     material_props = material_strings_.at(material_key);
    return material_props;
  }

  int
  GetBlockIdFromMaterial(const std::string& material_key) const
  {
    int block_id = -1;
    // std::map<int, BlockProperties>::iterator it;
    // for(it = macroscale_blocks_.begin(); it != macroscale_blocks_.end();
    // it++){
    for (auto const& entry : macroscale_blocks_) {
      BlockProperties const& block_props = entry.second;  // it->second;
      if (material_key == block_props.material_key_) {
        block_id = entry.first;  // it->first; //Found the block id with the
                                 // material string
        break;
      }
    }
    return block_id;
  }

  std::map<int, std::string>
  GetMicroscaleMaterialParameters() const
  {
    std::map<int, std::string> block_id_to_microscale_props_map;
    for (auto const& entry : microscale_blocks_) {
      int                    block_id            = entry.first;
      BlockProperties const& block_props         = entry.second;
      std::string const&     material_props      = material_strings_.at(block_props.material_key_);
      block_id_to_microscale_props_map[block_id] = material_props;
    }
    return block_id_to_microscale_props_map;
  }

  std::string
  GetMicroscaleBoundaryConditionStrategy() const
  {
    return microscale_boundary_condition_strategy_;
  }

  std::vector<std::string> const&
  GetBoundaryConditionStrings() const
  {
    return boundary_condition_strings_;
  }

  std::string
  GetOutputFieldString() const
  {
    if (output_field_string_.empty()) {
      std::string msg =
          "\n**** Error in Parser::GetOutputFieldString(), output fields not "
          "found (possible input deck error?).";
      throw std::logic_error(msg);
    }
    return output_field_string_;
  }

  /// \brief Set that Kokkos is used for the simulation
  ///
  /// \note The function will throw an exception when an environment
  /// has been set previously
  void
  SetToUseKokkos()
  {
    if ((env_set_) && (!use_kokkos_))
      throw std::runtime_error(" Conflicting Environment Variable ");
    use_kokkos_ = true;
    env_set_    = true;
  }

  /// \brief Set that Tpetra is used for the simulation
  ///
  /// \note The function will throw an exception when an environment
  /// has been set previously
  void
  SetToUseTpetra()
  {
    if ((env_set_) && (!use_tpetra_))
      throw std::runtime_error(" Conflicting Environment Variable ");
    use_tpetra_ = true;
    env_set_    = true;
  }

  /// \brief Set that UQ is used for the simulation
  ///
  void
  SetToUseUQ()
  {
    use_uq_  = true;
  }

  /// \brief Set that VT is used for the simulation
  ///
  /// \note The function will throw an exception when an environment
  /// has been set previously
  void
  SetToUseVT()
  {
    if ((env_set_) && (!use_vt_))
      throw std::runtime_error(" Conflicting Environment Variable ");
    use_vt_  = true;
    env_set_ = true;
  }

  bool
  UseKokkos() const
  {
    return use_kokkos_;
  }

  /// \brief Indicate whether Tpetra is used for the simulation
  bool
  UseTpetra() const
  {
    return use_tpetra_;
  }

  /// \brief Indicate whether UQ is used for the simulation
  bool
  UseUQ() const
  {
    return use_uq_;
  }

  /// \brief Indicate whether VT is used
  bool
  UseVT() const
  {
    return use_vt_;
  }

  /// \brief Set the rank index
  /// \param[in] myrank Rank index
  void
  SetRankID(int myrank)
  {
    my_rank_ = myrank;
  }

  /// \brief Get the rank index
  int
  GetRankID() const
  {
    return my_rank_;
  }

  /// \brief Set the number of ranks
  /// \param[in] nrank Number of MPI ranks
  void
  SetNumRanks(int nrank)
  {
    num_ranks_ = nrank;
  }

  /// \brief Get the number of ranks
  int
  GetNumRanks() const
  {
    return num_ranks_;
  }

  /// \brief Set the string for the input file
  /// \param[in] name Input filename
  void
  SetInputFilename(const std::string& name)
  {
    file_name_ = name;
  }

#ifdef NIMBLE_HAVE_TRILINOS
  void
  ResetTpetraScope(Tpetra::ScopeGuard* ptr_scope)
  {
    tpetra_scope_.reset(ptr_scope);
  }
#endif

#ifdef NIMBLE_HAVE_BVH
  int
  ContactDicing() const noexcept
  {
    return contact_dicing_;
  }
#endif

#ifdef NIMBLE_HAVE_UQ
  std::string const&
  UqModelString() const
  {
    return uq_model_string_;
  }
  std::map<std::string, std::string> const&
  UqParamsStrings() const
  {
    return uq_parameters_strings_;
  }

  bool
  HasUq() const
  {
    if (uq_model_string_.size() == 0 || uq_parameters_strings_.size() == 0)
      return false;
    else
      return true;
  }
#endif

 protected:
  virtual void
  ParseKeyValue(const std::string& k, const std::string& v);

  void
  ReadFile();

  std::string                        genesis_file_name_;
  std::string                        rve_genesis_file_name_;
  std::string                        exodus_file_name_;
  bool                               use_two_level_mesh_decomposition_;
  bool                               write_timing_data_file_;
  double                             nonlinear_solver_relative_tolerance_;
  int                                nonlinear_solver_max_iterations_;
  std::string                        time_integration_scheme_;
  double                             initial_time_{0.0};
  double                             final_time_{0.0};
  int                                num_load_steps_;
  int                                output_frequency_;
  std::string                        contact_string_;
  std::string                        contact_backend_string_;
  bool                               visualize_contact_entities_;
  bool                               visualize_contact_bounding_boxes_;
  std::string                        contact_visualization_file_name_;
  std::vector<int>                   microscale_output_element_ids_;
  std::map<std::string, std::string> material_strings_;
  std::map<int, BlockProperties>     macroscale_blocks_;
  std::map<int, BlockProperties>     microscale_blocks_;
  std::string                        microscale_boundary_condition_strategy_;
  std::vector<std::string>           boundary_condition_strings_;
  std::string                        output_field_string_;

#ifdef NIMBLE_HAVE_BVH
  int contact_dicing_ = 1;
#endif

#ifdef NIMBLE_HAVE_UQ
  std::map<std::string, std::string> uq_parameters_strings_;
  std::string                        uq_model_string_;
#endif

  //--- Variables related to environment

  /// \brief Boolean indicating whether the environment has been set
  bool env_set_ = false;

  /// \brief Filename for input file to read
  std::string file_name_ = std::string("none");

  /// \brief Boolean setting the usage of VT
  bool use_vt_ = false;

  /// \brief Boolean setting the usage of Kokkos
  bool use_kokkos_ = false;

  /// \brief Boolean setting the usage of Tpetra
  bool use_tpetra_ = false;

  /// \brief Boolean setting the usage of UQ
  bool use_uq_ = false;

  /// \brief Rank ID when running with MPI
  int my_rank_ = 0;

  /// \brief Number of ranks when running with MPI
  int num_ranks_ = 1;

#ifdef NIMBLE_HAVE_TRILINOS
  std::unique_ptr<Tpetra::ScopeGuard> tpetra_scope_ = nullptr;
#endif
};
}  // namespace nimble

#endif  // NIMBLE_PARSER_H
