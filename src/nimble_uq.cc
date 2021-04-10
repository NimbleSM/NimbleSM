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
      seed_(0),
      analyze_data_(false),
      write_text_data(false)
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
// parses uq model: quadature/samples {N} {N_exact} seed {SEED}
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
  } else {
    //    throw std::logic_error("\nUQ: unsupported !!\n");
    nsamples_ = std::stoi(items.front());
    items.pop_front();
    nexact_samples_ = std::stoi(items.front());
    items.pop_front();
    // sampling seed options
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

  if (items.size() > 0) {
    write_text_data = true;
    printf("WARNING: UQ will write text files\n");
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
    std::map<int, nimble_uq::Block*>&                   blocks)
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
  if (block_id_present) { blocks[block_id]->SetUqParameters(indices); }
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
  } else {
    SampleParameters();  // NEED TO DO THIS ONLY IF THEY HAVENT BEEN READ FROM FILE
    CenterParameters();
    WriteSamples();
  }
  Allocate();
  InitializeAnalyses();
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
  int n = GetNumSamples();
  int offnominal_displacement_ids[n],  // tmp?
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
UqModel::UpdateVelocity(double dt)
{
  if (!initialized_) { return; }
  int n = GetNumSamples();
  for (int s = 0; s < n; ++s) {
    double* f = offnominal_forces_[s];
    double* v = offnominal_velocities_[s];
    for (int i = 0; i < nunknowns_; ++i) {
      double m = mass_[i / 3];
      double a = (1.0 / m) * f[i];
      v[i] += dt * a;
    }
  }
}
//===========================================================================
void
UqModel::UpdateDisplacement(double dt)
{
  if (!initialized_) { return; }
  int n = GetNumSamples();
  for (int s = 0; s < n; ++s) {
    double* u = offnominal_displacements_[s];
    double* v = offnominal_velocities_[s];
    for (int i = 0; i < nunknowns_; ++i) { u[i] += dt * v[i]; }
  }
}
//===========================================================================
void
UqModel::Prep()
{
  if (!initialized_) { return; }
  for (int s = 0; s < nsamples_; ++s) {
    double* f = offnominal_forces_[s];
    for (int i = 0; i < nunknowns_; ++i) { f[i] = 0.0; }
  }
}
//===========================================================================
void
UqModel::ApplyClosure(const double* const nominal_internal_force)
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
      f[i] = nominal_internal_force[i] * interpolation_coefficients_[s][0];  // 0.0;
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
void
UqModel::InitializeAnalyses()
{
  if (!analyze_data_) { return; }
  if (!initialized_) { return; }
  int ns = GetNumSamples();
  max_V_stress_.reserve(ns);
  max_T_stress_.reserve(ns);
  std::fill(max_V_stress_.begin(), max_V_stress_.end(), 0.0);  // Initialize to 0.0
  std::fill(max_T_stress_.begin(), max_T_stress_.end(), 0.0);  // Initialize to 0.0
  FILE* fp;
  fp = fopen("QoIs.dat", "w");
  fclose(fp);
}
//===========================================================================
void
UqModel::WriteQoIs()
{  // currently not called
  if (!analyze_data_) { return; }
  if (!initialized_) { return; }
  int   ns              = GetNumSamples();
  int   nexacts         = GetNumExactSamples();
  int   napprox_samples = ns - nexacts;
  FILE* fp;
  fp = fopen("QoIs.dat", "w+");
  for (int nsmpl = 0; nsmpl < ns; nsmpl++) {
    for (int np = 0; np < nparameters_; np++) { fprintf(fp, "%e ", parameter_samples_[nsmpl][np]); }
    fprintf(fp, "%e %e %e %e \n", max_V_stress_[nsmpl], max_T_stress_[nsmpl], 0.0, 0.0);
  }
  fprintf(fp, "\n\n");
  fclose(fp);
}
//===========================================================================
void
UqModel::PerformAnalyses(
    const double* const reference_coordinates,
    int                 num_elem,
    const int* const    elem_conn,
    int                 block_id,
    Block&              block)
{  // called per block
  if (!analyze_data_) { return; }
  if (!initialized_) { return; }
  int dim                 = block.GetElementPointer()->Dim();
  int num_node_per_elem   = block.GetElementPointer()->NumNodesPerElement();
  int num_int_pt_per_elem = block.GetElementPointer()->NumIntegrationPointsPerElement();

  int vector_size      = LengthToInt(VECTOR, dim);
  int full_tensor_size = LengthToInt(FULL_TENSOR, dim);
  int sym_tensor_size  = LengthToInt(SYMMETRIC_TENSOR, dim);

  double ref_coord[vector_size * num_node_per_elem];
  double cur_coord[vector_size * num_node_per_elem];
  double def_grad[full_tensor_size * num_int_pt_per_elem];
  double cauchy_stress[sym_tensor_size * num_int_pt_per_elem];

  double vol_ave_stress[sym_tensor_size];
  double volume;

  // Find the indices of bulk_modulus and shear_modulus for this block from the parameters list
  // The parameters names span all blocks, so the names "bulk_modulus" could repeat across blocks
  int                                first_param = block_first_param_index_[block_id];
  int                                last_param  = block_last_param_index_[block_id];
  std::vector<std::string>::iterator first, last, found;
  first = names_.begin() + first_param;
  last  = names_.begin() + last_param + 1;

  found              = std::find(first, last, "bulk_modulus");
  int bulk_mod_index = std::distance(names_.begin(), found);

  found               = std::find(first, last, "shear_modulus");
  int shear_mod_index = std::distance(names_.begin(), found);

  int ns      = GetNumSamples();
  int nexacts = GetNumExactSamples();
  for (int sample_to_analyse = nexacts; sample_to_analyse < ns;
       sample_to_analyse++) {  // TESTING ON LAST EXACT SAMPLE TRAJECTORY
    double*             displacement    = offnominal_displacements_[sample_to_analyse];
    std::vector<double> material_params = parameter_samples_[sample_to_analyse];

    for (int elem = 0; elem < num_elem; elem++) {
      for (int node = 0; node < num_node_per_elem; node++) {
        int node_id = elem_conn[elem * num_node_per_elem + node];
        for (int i = 0; i < vector_size; i++) {
          ref_coord[node * vector_size + i] = reference_coordinates[vector_size * node_id + i];
          cur_coord[node * vector_size + i] =
              reference_coordinates[vector_size * node_id + i] + displacement[vector_size * node_id + i];
        }
      }

      block.GetElementPointer()->ComputeDeformationGradients(ref_coord, cur_coord, def_grad);

      block.GetMaterialPointer()->GetOffNominalStress(
          material_params[bulk_mod_index],
          material_params[shear_mod_index],
          num_int_pt_per_elem,
          def_grad,
          cauchy_stress);

      block.GetElementPointer()->ComputeVolumeAverage(
          cur_coord, sym_tensor_size, cauchy_stress, volume, vol_ave_stress);

      double max_elem_stress = 0.0;

      double evec1[3] = {0.0, 0.0, 0.0};  // Eigen vectors of cauchy stress
      double evec2[3] = {0.0, 0.0, 0.0};
      double evec3[3] = {0.0, 0.0, 0.0};

      double eval[3] = {0.0, 0.0, 0.0};  // Eigen values of cauchy stress

      Eigen_Sym33_NonUnit(vol_ave_stress, eval[0], eval[1], eval[2], evec1, evec2, evec3);

      double eigen_stress = 0.0;
      if (block_id == 1) {  // Von Mises stress
        eigen_stress = (eval[0] - eval[1]) * (eval[0] - eval[1]);
        eigen_stress += (eval[1] - eval[2]) * (eval[1] - eval[2]);
        eigen_stress += (eval[2] - eval[0]) * (eval[2] - eval[0]);

        eigen_stress = sqrt(eigen_stress);
      } else if (block_id == 2) {  // Tensile stress
        eigen_stress = *std::max_element(eval, eval + 3);
      } else {
        throw std::logic_error("\nUQ: Error unknown block id " + std::to_string(block_id) + " in UQ analyses !!!\n");
      }

      // Set the max stress in this elem
      max_elem_stress = std::max(eigen_stress, max_elem_stress);

      // Set the runnimg max over all elements for this uq sample
      if (block_id == 1) {  // Von Mises stress
        max_V_stress_[sample_to_analyse] = std::max(max_elem_stress, max_V_stress_[sample_to_analyse]);
      } else if (block_id == 2) {  // Tensile stress
        max_T_stress_[sample_to_analyse] = std::max(max_elem_stress, max_T_stress_[sample_to_analyse]);
      }

    }  // End loop over elem

  }  // End loop over samples_to_analyse

  // HKprintf("Max stress in block %d this step is %e\n", block_id, (block_id==1) ? max_V_stress_[sample_to_analyse] :
  // max_T_stress_[sample_to_analyse] );
}
//===========================================================================
void
UqModel::Write(int step)
{
  if (!initialized_) { return; }
  if (!write_text_data) { return; }
  int vector_size                   = LengthToInt(VECTOR, ndims_);
  int full_tensor_size              = LengthToInt(FULL_TENSOR, ndims_);
  int sym_tensor_size               = LengthToInt(SYMMETRIC_TENSOR, ndims_);
  int reference_coordinate_field_id = data_->GetFieldId("reference_coordinate");

  FILE*   fp;
  double* reference_coordinates = data_->GetNodeData(reference_coordinate_field_id);

  int ns      = GetNumSamples();
  int nexacts = GetNumExactSamples();

  std::map<int, nimble::Block>&          blocks = data_->GetBlocks();
  std::map<int, nimble::Block>::iterator block_it;

  std::string fname;
  for (int sample_to_analyse = nexacts; sample_to_analyse < ns; sample_to_analyse++) {
    int                 idx             = sample_to_analyse - nexacts;
    double*             displacement    = offnominal_displacements_[sample_to_analyse];
    std::vector<double> material_params = parameter_samples_[sample_to_analyse];
    fname                               = "displacement" + std::to_string(idx) + "_" + std::to_string(step) + ".dat";
    fp                                  = fopen(fname.c_str(), "w");
    for (int i = 0; i < nnodes_; i++) {
      int    ii = vector_size * i;
      double ux = displacement[ii];
      double uy = displacement[ii + 1];
      fprintf(fp, "%e %e\n", ux, uy);  // 2D
    }
    fclose(fp);
    fname = "stress" + std::to_string(idx) + "_" + std::to_string(step) + ".dat";
    fp    = fopen(fname.c_str(), "w");
    for (block_it = blocks.begin(); block_it != blocks.end(); block_it++) {
      int            block_id            = block_it->first;
      nimble::Block& block               = block_it->second;
      int            num_elem            = mesh_->GetNumElementsInBlock(block_id);
      int            num_node_per_elem   = block.GetElementPointer()->NumNodesPerElement();
      int            num_int_pt_per_elem = block.GetElementPointer()->NumIntegrationPointsPerElement();
      // scratch
      double ref_coord[vector_size * num_node_per_elem];
      double cur_coord[vector_size * num_node_per_elem];
      double def_grad[full_tensor_size * num_int_pt_per_elem];
      double cauchy_stress[sym_tensor_size * num_int_pt_per_elem];
      double vol_ave_stress[sym_tensor_size];
      double volume;
      // material parameters
      int                                first_param = block_first_param_index_[block_id];
      int                                last_param  = block_last_param_index_[block_id];
      std::vector<std::string>::iterator first, last, found;
      first = names_.begin() + first_param;
      last  = names_.begin() + last_param + 1;
      // Find the indices of bulk_modulus and shear_modulus for this block from the parameters list
      // The parameters names span all blocks, so the names "bulk_modulus" could repeat across blocks
      found               = std::find(first, last, "bulk_modulus");
      int bulk_mod_index  = std::distance(names_.begin(), found);
      found               = std::find(first, last, "shear_modulus");
      int shear_mod_index = std::distance(names_.begin(), found);

      int const* elem_conn = mesh_->GetConnectivity(block_id);
      for (int elem = 0; elem < num_elem; elem++) {
        for (int node = 0; node < num_node_per_elem; node++) {
          int node_id = elem_conn[elem * num_node_per_elem + node];
          for (int i = 0; i < vector_size; i++) {
            ref_coord[node * vector_size + i] = reference_coordinates[vector_size * node_id + i];
            cur_coord[node * vector_size + i] =
                reference_coordinates[vector_size * node_id + i] + displacement[vector_size * node_id + i];
          }
        }

        block.GetElementPointer()->ComputeDeformationGradients(ref_coord, cur_coord, def_grad);

        block.GetMaterialPointer()->GetOffNominalStress(
            material_params[bulk_mod_index],
            material_params[shear_mod_index],
            num_int_pt_per_elem,
            def_grad,
            cauchy_stress);

        block.GetElementPointer()->ComputeVolumeAverage(
            cur_coord, sym_tensor_size, cauchy_stress, volume, vol_ave_stress);
        double sxx = vol_ave_stress[K_S_XX];
        double syy = vol_ave_stress[K_S_YY];
        double sxy = vol_ave_stress[K_S_XY];
        fprintf(fp, "%e %e %e\n", sxx, syy, sxy);  // 2D

      }  // End loop over elem
    }    // End loop over blocks
    fclose(fp);
  }  // End loop over samples_to_analyse
}
//===========================================================================
}  // namespace nimble
