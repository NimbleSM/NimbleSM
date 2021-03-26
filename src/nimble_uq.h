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

#ifndef NIMBLE_UQ_H
#define NIMBLE_UQ_H

#include <memory>
#include <string>
#include <vector>
#include <map>

namespace nimble {
  class GenesisMesh;;
  class ModelData;
  class Block;
  class MaterialFactory;

//========================================================================
  class UqModel { // model
//========================================================================

  public:

    UqModel(int ndim, int nnodes, int nblocks=1);
    virtual ~UqModel() = default;

    void ParseConfiguration(
      std::string const & line);

    void ParseBlockInput(
      std::string const & line, 
      const int & block_id,
      std::string const & nominal_params_string,
      MaterialFactory& material_factory,
      bool block_id_present,
      std::map<int, Block> & blocks);

    // pre time integration
    void Initialize(GenesisMesh & mesh, ModelData & data);
    // post time integration
    void Finalize();
    // where all memory allocation is done
    void Allocate();
    // pre ComputeInternalForce
    void Setup();
    void Prep();
    // trajectories
    void UpdateVelocity(double dt); 
    void UpdateDisplacement(double dt); 
    void ApplyClosure(const double * const nominal_internal_force);//Do the exact_forces --> approx_forces map
//  void ApplyBoundaryConditions();

    // accessors
    int GetNumNodes()
      { return nnodes_;}
    int GetNumDims()
      { return ndims_;}
    std::vector<std::string> const & GetParameterNames() const 
      { return names_; }
    std::vector<double> const & GetParameters(const int & sample_id ) const 
       { return parameter_samples_[sample_id] ; }
    int GetNumSamples() const 
      { return nsamples_; }
    int GetNumParameters() const 
      { return nominal_parameter_values_.size(); }
    int GetNumExactSamples() const
      { return nexact_samples_; }
    // - pass throughs
    bool Initialized() const 
      { return initialized_; }
    std::vector<double*> & Displacements() 
      { return offnominal_displacements_ ; }
    std::vector<double*> & Velocities() 
       { return offnominal_velocities_ ; }
    std::vector<double*> & Forces() 
       { return offnominal_forces_ ; }

    // initialization 
    // - sample given distributions
    void SampleParameters();
    // - read samples
    void ReadParameters(std::string filename){}; // NOT IMPLEMENTED YET
    // - center parameter values on nominal
    void CenterParameters();
    // - input
    void ReadSamples();
    // - scale input
    void ScaleParameters();
    // - output
    void WriteSamples();
    // accessors
    std::vector<double> & GetNominalParameters() 
      { return  nominal_parameter_values_;} // NOTE 1 nominal trajectory
    std::vector<std::vector<double>> & GetParameterSamples()
      { return  parameter_samples_; }
    std::vector<double> & GetRanges() 
      { return  ranges_; }
    std::vector<std::vector<double>> GetParameterDifferences()
      { return  parameter_differences_; }

    //Routines for in-situ/in-line analyses over UQ trajectories
    void InitializeAnalyses();

    //Routine to perform in-situ processing over sample trajectories
    void PerformAnalyses(const double * const reference_coordinates, 
                         int num_elem, const int * const elem_conn,
                         int block_id, Block & block);
    // write current QoIs
    void WriteQoIs();

    // write data
    void Write(int step);

  private:
    int ndims_;  // from mesh
    int nnodes_;  // from mesh
    int nblocks_; // from mesh
    int nunknowns_;
    // parameter data
    std::string sampling_strategy_; 
    bool samples_from_file_;
    int seed_;
    int nsamples_; // number of _requested_ samples not necessarily all
    int nexact_samples_; //Subset of nsamples_ for which exact trajectories are computed
    int nparameters_;
    GenesisMesh * mesh_;
    ModelData * data_;
    const double* mass_ ;
    std::string samples_fname_;
    std::vector<std::string> names_;
    std::vector<std::string> distributions_; // independent distributions 
    std::vector<double> ranges_; // width 
    std::vector<double> nominal_parameter_values_; // NOTE 1 nominal trajectory
    std::vector<std::vector<double>> parameter_samples_;     // ns X np
    std::vector<std::vector<double>> parameter_differences_; // ns X np
    std::map<int, int> block_first_param_index_;
    std::map<int, int> block_last_param_index_;

    std::vector<std::vector<double>> interpolation_coefficients_; //napprox_smpls X np

    std::vector<double*> offnominal_displacements_   ;
    std::vector<double*> offnominal_velocities_      ;
    std::vector<double*> offnominal_forces_          ;

    //FIX THIS: Currently hardcoding for a specific analyses; Generalize this
    std::vector<double> max_V_stress_, max_T_stress_;

    bool initialized_;
    bool analyze_data_;
    bool write_text_data;
  };

} // namespace nimble

#endif // NIMBLE_UQ_H
