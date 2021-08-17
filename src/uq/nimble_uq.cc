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

#include "nimble_uq.h"

#include <algorithm>
#include <fstream>
#include <list>
#include <random>
#include <sstream>
#include <string>

#include "nimble_block.h"
#include "nimble_data_manager.h"
#include "nimble_data_utils.h"
#include "nimble_material_factory.h"
#include "nimble_uq_block.h"
#include "nimble_utils.h"

namespace nimble {
//===========================================================================
// NOTE could move this to nimble_utils
std::list<std::string>
split(std::string line)
{
  auto iss   = std::istringstream{line};
  auto str   = std::string{};
  auto items = std::list<std::string>();
  while (iss >> str) { items.push_back(str); }
  return items;
};

//===========================================================================
UqModel::UqModel(int ndims, const GenesisMesh* mesh, ModelData* data)
    : ndims_(ndims),
      nnodes_(mesh->GetNumNodes()),
      nblocks_(mesh->GetNumBlocks()),
      nsamples_(0),
      nparameters_(0),
      mesh_(mesh),
      data_(data),
      nexact_samples_(0)
{
  nunknowns_ = ndims_ * nnodes_;
  names_.clear();
  distributions_.clear();
  ranges_.clear();
  nominal_parameter_values_.clear();
  parameter_samples_.clear();
  interpolation_coefficients_.clear();
  initialized_       = false;
}
//===========================================================================
// parses 
//        uq model: file {FNAME}
void
UqModel::ParseConfiguration(std::string const& line)
{
  std::list<std::string> items = split(line);
  std::string sampling = items.front();
  items.pop_front();
  if (sampling == "file") {
    samples_fname_ = items.front();
    items.pop_front();
  } else { 
    std::cerr << " !! Error when parsing UQ parameters : !! \n";
  }
//std::cout << " reading samples from " << samples_fname_ << "\n";
}
//===========================================================================
// parses: 
// uq parameters: material_1 {PNAME1A} {STDDEV1A}
// uq parameters: material_2 {PNAME2A} {STDDEV1A} {PNAME2B} {STDDEV1B}
// uq parameters: material_2 {PNAME2A} {STDDEV1A} {PNAME2B} {STDDEV1B}
void
UqModel::ParseBlockInput(
    std::string const&                                  line,
    const int&                                          block_id,
    std::string const&                                  nominal_params_string,
    const std::shared_ptr<nimble::MaterialFactoryBase>& material_factory,
    bool                                                block_id_present,
    std::map<int, std::shared_ptr<nimble::Block>>&      blocks)
{
  std::cout << "START ParseBlockInput\n" << std::flush;
  if (line == "none")
    throw std::logic_error(
        "\nError: NIMBLE_DO_UQ is defined but no uq parameters are specified for material corresponding to block id " +
        std::to_string(block_id) + "\n");
  int                    start_idx = GetNumParameters();
  std::list<std::string> items     = split(line);
  while (items.size() > 0) {
      std::string pname = items.front();
      items.pop_front();
      double range = std::stod(items.front());
      items.pop_front();
      names_.push_back(pname);
      ranges_.push_back(range);
      nparameters_++;
  }

  std::map<std::string, double> nominal_material_params =
      material_factory->parse_material_params_string(nominal_params_string);
  std::map<std::string, int> indices;
  for (int idx = start_idx; idx < nparameters_; idx++) {
    std::string pname = names_[idx];
    indices[pname]    = idx;
    if (pname == "density") {
      nominal_parameter_values_.push_back(nominal_material_params[pname]);  // HK( block.GetDensity() );
    } else if (pname == "bulk_modulus") {
      nominal_parameter_values_.push_back(nominal_material_params[pname]);  // HK( block.GetBulkModulus() ); <<<<?????
    } else if (pname == "shear_modulus") {
      nominal_parameter_values_.push_back(nominal_material_params[pname]);  // HK( block.GetShearModulus() );
    } else {
      throw std::logic_error("\nError: uq parameter " + pname + " is unrecognised or not currently supported\n");
    }
  }
  block_first_param_index_[block_id] = start_idx;
  block_last_param_index_[block_id]  = nparameters_ - 1;
  // let block know if block id is present on this rank
  if (block_id_present) {
    auto uq_block_ptr = dynamic_cast<nimble_uq::Block*>(blocks[block_id].get());
    uq_block_ptr->SetUqParameters(indices);
  }
  std::cout << "FINISH ParseBlockInput\n" << std::flush;
}
//===========================================================================
void
UqModel::ReadSamples()
{
//std::cout << " start samples from " << samples_fname_ << "\n";
  std::ifstream          f(samples_fname_);
  std::string            line;
  std::list<std::string> items;

  // read 1st line
  int linenumber = 0;
  getline(f, line); linenumber++;
  items.clear();
  items = split(line);
  nexact_samples_ = std::stoi(items.front()); // in addition to the nominal
  items.pop_front();
  napprox_samples_ = std::stoi(items.front());
  items.pop_front();
  nsamples_ = napprox_samples_ + nexact_samples_;

  parameter_samples_.resize(nsamples_);
  interpolation_coefficients_.resize(napprox_samples_);

  //- (A) -the exact samples (no interp coefficients for these)
  //--{exact_smpl_1_P1} {exact_smpl_1_P2}
  //--{exact_smpl_2_P1} {exact_smpl_2_P2}
  for (int nsmpl = 0; nsmpl < nexact_samples_; nsmpl++) {
    getline(f, line); linenumber++;
    items.clear();
    items = split(line);
    if (items.size() != nparameters_) {
      throw std::logic_error("\nError: UQ file " + samples_fname_ + " has not enough exact sample parameters on line " + std::to_string(linenumber) + " expecting " + std::to_string(nparameters_) + " found " + std::to_string(items.size()) + "\n");
    }
    for (int n = 0; n < nparameters_; n++) {
      std::string item = items.front();
      items.pop_front();
      parameter_samples_[nsmpl].push_back(std::stod(item));  // exact samples come before the approx samples
    }
  }
  std::cout << " read " << nexact_samples_ << " exact samples\n";
  //-- (B) the approximate samples, and interpolation coefficients
  //--{approx_smpl_1_P1} {approx_smpl_1_P2} : {interp_coeff_1} {interp_coeff_2}
  //--{approx_smpl_2_P1} {approx_smpl_2_P2} : {interp_coeff_1} {interp_coeff_2}
  int ncoeff = nexact_samples_ + 1; // including nominal
  for (int nsmpl = 0; nsmpl < napprox_samples_; nsmpl++) {
    getline(f, line); linenumber++;
    size_t colon_pos = line.find(":");
    if (colon_pos == std::string::npos) {
      throw std::logic_error( "\nError: UQ file " + samples_fname_ + " has wrong format (no colon) on line " + std::to_string(linenumber) + "\n");
    }
    std::string params_string = line.substr(0, colon_pos);
    std::string coeffs_string = line.substr(colon_pos + 1, line.size());
    items.clear();
    items = split(params_string);
    if (items.size() != nparameters_) {
      throw std::logic_error( "\nError: UQ file " + samples_fname_ + " has " + std::to_string(items.size()) +" parameters (preceding :) on line " + std::to_string(linenumber) + " expecting " + std::to_string(nparameters_) + "\n");
    }
    for (int n = 0; n < nparameters_; n++) {
      std::string item = items.front();
      items.pop_front();
      parameter_samples_[nexact_samples_ + nsmpl].push_back(std::stod(item));  // approx samples come after the exact samples
    }
    items.clear();
    items = split(coeffs_string);
    if (items.size() != ncoeff) {
      throw std::logic_error( "\nError: UQ file " + samples_fname_ + " has " + std::to_string(items.size()) + " coefficients on line " + std::to_string(linenumber) + " expecting " + std::to_string(ncoeff) + "\n");
    }

    for (int n = 0; n < ncoeff; n++) {
      std::string item = items.front();
      items.pop_front();
      interpolation_coefficients_[nsmpl].push_back(std::stod(item));
    }
  }
  std::cout << " read " << napprox_samples_ << " approximate samples\n";
}
//===========================================================================
void
UqModel::ScaleParameters()
{
  //--Parameter values are in  -1.0:1.0 range
  //--Scaling them to physical units 
  //             -1.0 ==> nominal - range
  //              1.0 ==> nominal + range
  for (int p = 0; p < nparameters_; p++) {
    double slope = ranges_[p];
    double shift = nominal_parameter_values_[p];
    for (int s = 0; s < nsamples_; s++) {
      parameter_samples_[p][s] *= slope;
      parameter_samples_[p][s] += shift;
    }
  }
//std::cout << "done with ScaleParameters\n" << std::flush;
}
//===========================================================================
void
UqModel::WriteSamples()
{
  std::ofstream f("parameter_samples.dat");
  f << "# ";
  for (int p = 0; p < nparameters_; p++) { f << names_[p] << " "; }
  f << "\n";
  for (int s = 0; s < nsamples_; s++) {
    for (int p = 0; p < nparameters_; p++) { f << parameter_samples_[s][p] << " "; }
    f << "\n";
  }
  f.close();
//std::cout << "done with WriteSamples\n" << std::flush;
}
//===========================================================================
void
UqModel::Initialize()
{
  ReadSamples();
  ScaleParameters();  
  WriteSamples();
  Allocate();
  initialized_ = true;
  std::cout << "done with Initialize\n" << std::flush;
}
//===========================================================================
void
UqModel::Allocate()
{
  nominal_force_ = data_->GetVectorNodeData("internal_force").data(); 
  for (int i = 0; i < nsamples_; i++) {
    std::string tag = std::to_string(i+1);
    data_->AllocateNodeData(nimble::VECTOR, "sample_displacement_" + tag, nnodes_);
    data_->AllocateNodeData(nimble::VECTOR, "sample_velocity_" + tag, nnodes_);
    data_->AllocateNodeData(nimble::VECTOR, "sample_force_" + tag, nnodes_);
  }
  // for time averages
  data_->AllocateNodeData(nimble::VECTOR, "sum_displacement", nnodes_);
  data_->AllocateNodeData(nimble::VECTOR, "sum2_displacement", nnodes_);
  data_->AllocateNodeData(nimble::VECTOR, "min_displacement", nnodes_);
  data_->AllocateNodeData(nimble::VECTOR, "max_displacement", nnodes_);
}
//===========================================================================
void
UqModel::Setup()
{
  std::cout << "start with Setup\n" << std::flush;
  if (!initialized_) { return; }
  sample_displacements_.clear();
  sample_velocities_.clear();
  sample_forces_.clear();
// NOTE split into exact and samples?
  for (int i = 0; i < nsamples_; i++) {
    std::string tag = std::to_string(i+1);
    sample_displacements_.push_back( data_->GetNodeData(data_->GetFieldId("sample_displacement_" + tag)));
    sample_velocities_.push_back(    data_->GetNodeData(data_->GetFieldId("sample_velocity_"     + tag)));
    sample_forces_.push_back(        data_->GetNodeData(data_->GetFieldId("sample_force_"        + tag)));
  }
  std::cout << "done with Setup\n" << std::flush;
}
//===========================================================================
void
UqModel::ApplyClosure()
{
  if (!initialized_) { return; }
  // Now loop over the approximate samples, and apply the interpolation
  for (int s = 0; s < napprox_samples_; s++) {
    double* f = sample_forces_[s + nexact_samples_]; // approx after exact
    // Loop over each node
    for (int i = 0; i < nunknowns_; ++i) {
      f[i] = nominal_force_[i] * interpolation_coefficients_[s][0];  
      for (int k = 0; k < nexact_samples_; k++) {
        double* f_exact = sample_forces_[k];
        f[i] += f_exact[i] * interpolation_coefficients_[s][k + 1];
      }
    }
  }
}
//===========================================================================
}  // namespace nimble
