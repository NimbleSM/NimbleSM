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
      nexact_samples_(0),
      seed_(0)
{
  nunknowns_ = ndims_ * nnodes_;
  names_.clear();
  distributions_.clear();
  ranges_.clear();
  nominal_parameter_values_.clear();
  parameter_samples_.clear();
  parameter_differences_.clear();
  interpolation_coefficients_.clear();
  samples_from_file_ = false;
  initialized_       = false;
}
//===========================================================================
// parses uq model: quadature/samples {N} linear seed {SEED}
//        uq model: quadature/samples {N} {N_exact} seed {SEED}
//        uq model: file {N} {N_exact} {FNAME}
void
UqModel::ParseConfiguration(std::string const& line)
{
  std::list<std::string> items = split(line);
  // required
  sampling_strategy_ = items.front();
  items.pop_front();
  //  printf(sampling_strategy_.c_str(),"\n");

  if (sampling_strategy_ == "file") {
    samples_fname_ = items.front();
    items.pop_front();
    samples_from_file_ = true;
  } else { // NOTE REMOVE
    //    throw std::logic_error("\nUQ: unsupported !!\n");
    nsamples_ = std::stoi(items.front());
    items.pop_front();
    //
    std::string key = items.front();
    if (key == "linear") {
      items.pop_front();
    } else {
      try {
        nexact_samples_ = std::stoi(items.front());
      } catch (const std::invalid_argument& e) {
        throw;
      } catch (...) {
        std::cerr << " !! Error when parsing UQ parameters : " << items.front() << " !! \n";
      }
      items.pop_front();
    }
    //
    // sampling seed options
    //
    if (items.size() > 0) {
      std::string key = items.front();
      items.pop_front();
      if (key == "seed") {
        seed_ = std::stoi(items.front());
        items.pop_front();
      } else {
        throw std::logic_error("\nUQ: Error unknown keyword !!!\n");
      }
    }
    if (seed_ < 1) {
      std::random_device seeder;
      seed_ = seeder();
    }
    if (nexact_samples_ > nsamples_) {
      throw std::logic_error("\nUQ: Error. Specified num_exact_samples to be greater than num_all_samples !!!\n");
    }
  }
}
//===========================================================================
// parses:
// uq parameters: material_1 {PNAME1A} gaussian {S1A}
// uq parameters: material_2 {PNAME2A} gaussian {S1A} {PNAME2B} gaussian {S1B}
// uq parameters: material_2 {PNAME2A} {S1A} {PNAME2B} {S1B}
void
UqModel::ParseBlockInput(
    std::string const&                                  line,
    const int&                                          block_id,
    std::string const&                                  nominal_params_string,
    const std::shared_ptr<nimble::MaterialFactoryBase>& material_factory,
    bool                                                block_id_present,
    std::map<int, std::shared_ptr<nimble::Block>>&      blocks)
{
  if (line == "none")
    throw std::logic_error(
        "\nError: NIMBLE_DO_UQ is defined but no uq parameters are specified for material corresponding to block id " +
        std::to_string(block_id) + "\n");
  int                    start_idx = GetNumParameters();
  std::list<std::string> items     = split(line);
  if (samples_from_file_) {
    while (items.size() > 0) {
      std::string pname = items.front();
      items.pop_front();
      double range = std::stod(items.front());
      items.pop_front();
      names_.push_back(pname);
      ranges_.push_back(range);
      nparameters_++;
    }
  } else {
    while (items.size() > 0) {
      std::string pname = items.front();
      items.pop_front();
      std::string pdist = items.front();
      items.pop_front();
      double range = std::stod(items.front());
      items.pop_front();
      names_.push_back(pname);
      distributions_.push_back(pdist);  // tmp
      ranges_.push_back(range);         // tmp?
      nparameters_++;
    }
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
      nominal_parameter_values_.push_back(nominal_material_params[pname]);  // HK( block.GetBulkModulus() );
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
}
//===========================================================================
void
UqModel::ReadSamples()
{
  std::ifstream          f(samples_fname_);
  std::string            line;
  std::list<std::string> items;

  // read 1st line
  getline(f, line);
  //  printf(line.c_str());
  items.clear();
  items = split(line);
  //  printf(items.front().c_str());
  nsamples_ = std::stoi(items.front());
  items.pop_front();
  //  printf(items.front().c_str());
  nexact_samples_ = std::stoi(items.front());
  items.pop_front();
  //  printf("samples %d trajectories %d\n",nsamples_,nexact_samples_);

  int ns      = GetNumSamples();
  int nexacts = GetNumExactSamples();

  int napprox_samples = ns - nexacts;

  parameter_samples_.resize(ns);
  interpolation_coefficients_.resize(napprox_samples);
  // HKfor(int n=0; n<ns; n++) {
  // HK  parameter_samples_[n].resize(nparameters_);
  // HK}

  //--The file first has the approximate samples, and interpolation coefficients
  //--{approx_smpl_1_P1} {approx_smpl_1_P2} : {interp_coeff_1} {interp_coeff_2}
  //--{approx_smpl_2_P1} {approx_smpl_2_P2} : {interp_coeff_1} {interp_coeff_2}
  for (int nsmpl = 0; nsmpl < napprox_samples; nsmpl++) {
    getline(f, line);

    size_t colon_pos = line.find(":");
    if (colon_pos == std::string::npos) {
      throw std::logic_error(
          "\nError: UQ samples configuration file " + samples_fname_ + " has wrong format (no colon) on line " +
          std::to_string(nsmpl) + "\n");
    }
    std::string params_string = line.substr(0, colon_pos);
    std::string coeffs_string = line.substr(colon_pos + 1, line.size());

    items.clear();
    items = split(params_string);

    if (items.size() != nparameters_) {
      throw std::logic_error(
          "\nError: UQ samples configuration file " + samples_fname_ + " has not enough params (preceding :) on line " +
          std::to_string(nsmpl) + " expecting " + std::to_string(nparameters_) + " found " +
          std::to_string(items.size()) + "\n");
    }

    for (int n = 0; n < nparameters_; n++) {
      std::string item = items.front();
      items.pop_front();
      parameter_samples_[nexacts + nsmpl].push_back(std::stod(item));  // approx samples come after the exact samples
    }

    items.clear();
    items = split(coeffs_string);

    if (items.size() != nexacts + 1) {
      throw std::logic_error(
          "\nError: UQ samples configuration file " + samples_fname_ + " has " + std::to_string(items.size()) +
          " coefficients on line " + std::to_string(nsmpl) + " expecting " + std::to_string(nexacts + 1) + "\n");
    }

    for (int n = 0; n < nexacts + 1; n++) {
      std::string item = items.front();
      items.pop_front();
      interpolation_coefficients_[nsmpl].push_back(std::stod(item));
    }
  }

  //--Next an empty line
  getline(f, line);

  //--Next all the exact samples (no interp coefficients for these)
  //--{exact_smpl_1_P1} {exact_smpl_1_P2}
  //--{exact_smpl_2_P1} {exact_smpl_2_P2}
  for (int nsmpl = 0; nsmpl < nexacts; nsmpl++) {
    getline(f, line);

    items.clear();
    items = split(line);

    if (items.size() != nparameters_) {
      throw std::logic_error(
          "\nError: UQ samples configuration file " + samples_fname_ + " has not enough exact sample params on line " +
          std::to_string(nsmpl + napprox_samples) + " expecting " + std::to_string(nparameters_) + " found " +
          std::to_string(items.size()) + "\n");
    }

    for (int n = 0; n < nparameters_; n++) {
      std::string item = items.front();
      items.pop_front();
      parameter_samples_[nsmpl].push_back(std::stod(item));  // exact samples come before the approx samples
    }
  }
}
//===========================================================================
void
UqModel::ScaleParameters()
{
  //--Parameter Values specified in file are in a -1.0 : 1.0 Range
  //--Scaling them to absolute range

  //--ASSUMPTION: Samples are from whatever distribution with
  //             -1.0 ==> nominal - range
  //              1.0 ==> nominal + range
  //              e.g., if distribution is Gaussian, range=3*sigma contains 99.7% of samples

  int ns = GetNumSamples();

  for (int np = 0; np < nparameters_; np++) {
    double slope = ranges_[np];
    double shift = nominal_parameter_values_[np];

    for (int n = 0; n < ns; n++) {
      parameter_samples_[n][np] *= slope;
      parameter_samples_[n][np] += shift;
    }
  }
}
//===========================================================================
void
UqModel::Initialize()
{
  if (samples_from_file_) {
    ReadSamples();
    ScaleParameters();  // PARMETERS IN FILE ARE in -1.0 : 1.0 RANGE
  } else { // NOTE REMOVE THIS
    SampleParameters();  // NEED TO DO THIS ONLY IF THEY HAVENT BEEN READ FROM FILE
    CenterParameters();
    WriteSamples();
  }
  Allocate();
  initialized_ = true;
}
//===========================================================================
void
UqModel::Finalize()
{
}

//===========================================================================
void
UqModel::Allocate()
{
  nominal_force_ = data_->GetVectorNodeData("internal_force").data(); // NOTE
  int n = GetNumSamples();
  int offnominal_displacement_ids[n],  // tmp? NOTE
      offnominal_internal_force_ids[n], offnominal_velocity_ids[n];
  for (int nuq = 0; nuq < n; nuq++) {
    offnominal_displacement_ids[nuq] =
        data_->AllocateNodeData(nimble::VECTOR, "off_nom_displacement_" + std::to_string(nuq), nnodes_);
    offnominal_velocity_ids[nuq] =
        data_->AllocateNodeData(nimble::VECTOR, "off_nom_velocity_" + std::to_string(nuq), nnodes_);
    offnominal_internal_force_ids[nuq] =
        data_->AllocateNodeData(nimble::VECTOR, "off_nom_internal_force_" + std::to_string(nuq), nnodes_);
  }
  // for time averages
  int sum_displacement_id  = data_->AllocateNodeData(nimble::VECTOR, "sum_displacement", nnodes_);
  int sum2_displacement_id = data_->AllocateNodeData(nimble::VECTOR, "sum2_displacement", nnodes_);
  int min_displacement_id  = data_->AllocateNodeData(nimble::VECTOR, "min_displacement", nnodes_);
  int max_displacement_id  = data_->AllocateNodeData(nimble::VECTOR, "max_displacement", nnodes_);
}
//===========================================================================
void
UqModel::Setup()
{
  if (!initialized_) { return; }
  mass_ = data_->GetNodeData(data_->GetFieldId("lumped_mass"));
  offnominal_displacements_.clear();
  offnominal_velocities_.clear();
  offnominal_forces_.clear();
  for (int nuq = 0; nuq < nsamples_; nuq++) {
    offnominal_displacements_.push_back(
        data_->GetNodeData(data_->GetFieldId("off_nom_displacement_" + std::to_string(nuq))));
    offnominal_velocities_.push_back(data_->GetNodeData(data_->GetFieldId("off_nom_velocity_" + std::to_string(nuq))));
    offnominal_forces_.push_back(
        data_->GetNodeData(data_->GetFieldId("off_nom_internal_force_" + std::to_string(nuq))));
  }
}
//===========================================================================
void
UqModel::ApplyClosure()
{
  if (!initialized_) { return; }
  // TO DO: THIS NEEDS TO BE FLESHED OUT
  int ns      = GetNumSamples();
  int nexacts = GetNumExactSamples();

  int napprox_samples = ns - nexacts;

  // Now loop over the approximate samples, and apply the interpolation
  for (int s = 0; s < napprox_samples; s++) {
    double* f = offnominal_forces_[s + nexacts];  // The approx samples come after the exact
    // Loop over each node
    for (int i = 0; i < nunknowns_; ++i) {
      f[i] = nominal_force_[i] * interpolation_coefficients_[s][0];  // 0.0;
      for (int k = 0; k < nexacts; k++) {
        double* f_exact = offnominal_forces_[k];
        f[i] += f_exact[i] * interpolation_coefficients_[s][k + 1];
      }
    }
  }
}
//===========================================================================
void
UqModel::CenterParameters()
{
  int ns = GetNumSamples();
  parameter_differences_.resize(ns);
  for (int s = 0; s < ns; s++) {
    for (int p = 0; p < nparameters_; p++) {
      double lambda0 = nominal_parameter_values_[p];
      double lambda  = parameter_samples_[s][p];
      parameter_differences_[s].push_back(lambda - lambda0);
    }
  }
}
//===========================================================================
void
UqModel::WriteSamples()
{
  int           ns = GetNumSamples();
  std::ofstream f("parameter_samples.dat");
  f << "# ";
  for (int p = 0; p < nparameters_; p++) { f << names_[p] << " "; }
  f << "\n";
  for (int s = 0; s < ns; s++) {
    for (int p = 0; p < nparameters_; p++) { f << parameter_samples_[s][p] << " "; }
    f << "\n";
  }
  f.close();
}
//===========================================================================
void
UqModel::SampleParameters()
{
  // sample each parameter independently
  std::default_random_engine r(seed_);
  parameter_samples_.resize(nsamples_);
  // Loop over each of the parameters
  for (int p = 0; p < nparameters_; p++) {
    double      sample                  = 0.0;
    std::string param_distribution_name = distributions_[p];
    if (param_distribution_name == "uniform") {
      double low = nominal_parameter_values_[p] - ranges_[p];
      double hi  = nominal_parameter_values_[p] + ranges_[p];
      std::cout << "UQ: sampling uniform distribution for " << names_[p] << " in range " << low << "--" << hi
                << " with " << nsamples_ << " samples\n";
      std::uniform_real_distribution<double> uniform_dist(low, hi);
      for (int s = 0; s < nsamples_; s++) {
        sample = uniform_dist(r);
        parameter_samples_[s].push_back(sample);
      }
    } else if (param_distribution_name == "gaussian") {
      double mean  = nominal_parameter_values_[p];
      double stdev = ranges_[p];
      std::cout << "UQ: sampling normal/gaussian distribution for " << names_[p] << " with mean: " << mean
                << " stdev:" << stdev << "\n";
      std::normal_distribution<double> normal_dist(mean, stdev);
      for (int s = 0; s < nsamples_; s++) {
        sample = normal_dist(r);
        parameter_samples_[s].push_back(sample);
      }
    } else {  // unknown distribution, all_samples=nominal_value
      std::cout << "UQ: unsupported distribution " << distributions_[p] << " specified for " << names_[p]
                << " setting Dirac delta\n";
      for (int s = 0; s < nsamples_; s++) { parameter_samples_[s].push_back(nominal_parameter_values_[p]); }
    }
  }
}
//===========================================================================
}  // namespace nimble
