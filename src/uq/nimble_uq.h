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

#include <map>
#include <memory>
#include <string>
#include <vector>

namespace nimble {
class GenesisMesh;
class ModelData;
class Block;
class MaterialFactory;
class MaterialFactoryBase;

class UqModel
{ 
  

 public:
  UqModel(int ndim, const GenesisMesh* mesh, ModelData* data);
  virtual ~UqModel() = default;

  void
  ParseConfiguration(std::string const& line);

  void
  ParseBlockInput(
      std::string const&                                  line,
      const int&                                          block_id,
      std::string const&                                  nominal_params_string,
      const std::shared_ptr<nimble::MaterialFactoryBase>& material_factory,
      bool                                                block_id_present,
      std::map<int, std::shared_ptr<nimble::Block>>&      blocks);

  // pre time integration
  void
  Initialize();
  // post time integration
  void
  Finalize();
  // where all memory allocation is done
  void
  Allocate();
  // pre ComputeInternalForce
  void
  Setup();
  // construct approx_forces from the exact samples
  void
  ApplyClosure();
  // initialization
  void
  ReadSamples();
  // - scale input
  void
  ScaleParameters();
  // - write scaled parameters
  void
  WriteSamples();
  // Routine to perform in-situ processing over sample trajectories
  void
  PerformAnalyses(
    const double*                   reference_coordinates,
    int                             num_elem,
    const int*                      elem_conn,
    int                             block_id,
    std::shared_ptr<nimble::Block>& block) {};
  // write data
  void
  Write(double time) {};

  // accessors
  int
  GetNumSamples() const
  {
    return nsamples_;
  }
  int
  GetNumParameters() const
  {
    return parameter_order_.size();
  }
  int
  GetNumExactSamples() const
  {
    return nexact_samples_;
  }
  bool
  Initialized() const
  {
    return initialized_;
  }
  std::vector<double*>&
  Displacements()
  {
    return displacements_;
  }
  std::vector<double*>&
  Velocities()
  {
    return velocities_;
  }
  std::vector<double*>&
  Forces()
  {
    return forces_;
  }
  std::map<std::string,double>
  Parameters(int block_id, int sample_i) 
  {
    return parameters_[block_id][sample_i];
  }

 private:
  int ndims_;    // from mesh
  int nnodes_;   // from mesh
  int nblocks_;  // from mesh
  int nunknowns_; // nnodes*ndof
  // parameter data
  int                              nsamples_;        // total additional 
  int                              nexact_samples_;  
  int                              napprox_samples_; 
  int                              nparameters_;
  const GenesisMesh*               mesh_;
  ModelData*                       data_;
  std::string                      samples_file_name_;
  std::map<int,std::map<std::string,double>> nominal_parameters_; 
  std::map<int,std::map<std::string,double>> parameter_uncertainties_; 
  std::vector<std::pair<int,std::string>> parameter_order_;
  std::vector<std::vector<double>> parameter_samples_;  // ns,np
  std::vector<std::vector<double>> interpolation_coefficients_;  // na,np
  // trajectory data
  double* nominal_force_;
  std::vector<double*> displacements_;
  std::vector<double*> velocities_;
  std::vector<double*> forces_;
  std::map<int,std::vector<std::map<std::string,double>>> parameters_;
  bool initialized_;
};

}  // namespace nimble

#endif  // NIMBLE_UQ_H
