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

#include "nimble_parser.h"
#include <fstream>
#include <cstring>
#include <string>
#include <sstream>
#include <cstdlib>
#include <algorithm>

#include <iostream>

namespace nimble {

  std::string IOFileName(std::string const & serial_name,
                         std::string const & extension,
                         std::string const & label,
                         int my_rank,
                         int num_ranks) {

    if (serial_name == "none") {
      return serial_name;
    }

    std::string file_name = serial_name;
    size_t pos = file_name.rfind(".g");
    if (pos != std::string::npos){
      file_name = file_name.substr(0, pos);
    }
    pos = file_name.rfind(".e");
    if (pos != std::string::npos){
      file_name = file_name.substr(0, pos);
    }
    if (label.size() > 0) {
      file_name += "." + label;
    }
    file_name += "." + extension;
    if (num_ranks > 1) {
      std::stringstream num_ranks_ss;
      num_ranks_ss << num_ranks;
      std::stringstream my_rank_ss;
      my_rank_ss << my_rank;
      int num_zero_padding = num_ranks_ss.str().size() - my_rank_ss.str().size();
      std::stringstream ss;
      ss << file_name << "." << num_ranks << ".";
      for (int i = 0 ; i < num_zero_padding ; i++) {
        ss << 0;
      }
      ss << my_rank;
      //file_name += ss.str(); // DJL not thread safe??
      ss >> file_name;
    }
    return file_name;
  }

  void IOFileNameThreadSafe(char const * const serial_file_name,
                            char const * const extension,
                            char const * const label,
                            int my_mpi_rank,
                            int mpi_num_ranks,
                            int my_thread_rank,
                            int num_threads,
                            char * const file_name) {

    // copy the input file name to the output file name
    strcpy(file_name, serial_file_name);

    // strip off .g or .e suffix
    size_t len = strlen(file_name);
    if (len > 2) {
      file_name[len-2] = '\0';
    }

    // add the label, if any, to the file_name
    size_t label_len = strlen(label);
    if (label_len > 0) {
      strcat(file_name, ".");
      strcat(file_name, label);
    }

    // add the extension to the file name
    size_t extension_len = strlen(extension);
    if (extension_len > 0) {
      strcat(file_name, ".");
      strcat(file_name, extension);
    }

    if (mpi_num_ranks > 1) {

      // total number of mpi ranks
      char mpi_num_ranks_str[MAX_C_STR_SIZE];
      sprintf(mpi_num_ranks_str, "%d", mpi_num_ranks);
      strcat(file_name, ".");
      strcat(file_name, mpi_num_ranks_str);

      // padded mpi rank number
      char my_mpi_rank_str[MAX_C_STR_SIZE];
      sprintf(my_mpi_rank_str, "%d", my_mpi_rank);
      size_t required_str_len = strlen(mpi_num_ranks_str);
      size_t current_str_len = strlen(my_mpi_rank_str);
      strcat(file_name, ".");
      size_t num_zero_padding = required_str_len - current_str_len;
      for (size_t i=0 ; i<num_zero_padding ; i++) {
        strcat(file_name, "0");
      }
      strcat(file_name, my_mpi_rank_str);
    }

    if (num_threads > 1) {

      // total number of threads
      char num_threads_str[MAX_C_STR_SIZE];
      sprintf(num_threads_str, "%d", num_threads);
      strcat(file_name, ".");
      strcat(file_name, num_threads_str);

      // padded thread number
      char my_thread_rank_str[MAX_C_STR_SIZE];
      sprintf(my_thread_rank_str, "%d", my_thread_rank);
      size_t required_str_len = strlen(num_threads_str);
      size_t current_str_len = strlen(my_thread_rank_str);
      strcat(file_name, ".");
      size_t num_zero_padding = required_str_len - current_str_len;
      for (size_t i=0 ; i<num_zero_padding ; i++) {
        strcat(file_name, "0");
      }
      strcat(file_name, my_thread_rank_str);
    }
  }

  BlockProperties::BlockProperties(std::string props) {
    size_t space_pos = props.find(" ");
    block_name_ = props.substr(0, space_pos);
    material_key_ = props.substr(space_pos+1, props.size());
    size_t underscore_pos = block_name_.rfind("_");
    std::string block_id_str = block_name_.substr(underscore_pos+1, block_name_.size());
    std::stringstream ss;
    ss << block_id_str;
    ss >> block_id_;
  }

  Parser::Parser() : file_name_("none"), genesis_file_name_("none"),  use_two_level_mesh_decomposition_(false), write_timing_data_file_(false), rve_genesis_file_name_("none"), exodus_file_name_("none"),
                     time_integration_scheme_("explicit"), nonlinear_solver_relative_tolerance_(1.0e-6), nonlinear_solver_max_iterations_(200), final_time_(1.0), num_load_steps_(0), output_frequency_(1), visualize_contact_entities_(false), visualize_contact_bounding_boxes_(false), contact_visualization_file_name_("none"), microscale_boundary_condition_strategy_("periodic bc") {
    material_strings_["rve"] = "none";
  }

  void Parser::Initialize(std::string file_name) {
    file_name_ = file_name;
    ReadFile();
  }

  void Parser::ReadFile() {

    int max_chars_per_line = 1000;
    std::string delimiter(" ");

    std::ifstream fin;
    fin.open(file_name_.c_str());
    if ( !fin.good() ) {
      std::string msg = "\n**** Error in Parser::ReadFile(), unable to read file " + file_name_ + "\n";
      throw std::logic_error(msg);
    }

    while ( !fin.eof() ) {

      char buf[max_chars_per_line];
      fin.getline(buf, max_chars_per_line);
      std::string line(buf);

      // Strip out comments (# character)
      size_t pound_pos = line.find("#");
      if ( pound_pos != std::string::npos ) {
        line = line.substr(0, pound_pos);
      }
      // Strip out leading and trailing whitespace
      if ( line.size() > 0 ) {
        size_t start_pos = line.find_first_not_of(" \t");
        size_t end_pos = line.find_last_not_of(" \t");
        line = line.substr(start_pos, end_pos - start_pos + 1);
      }

      if ( line.size() > 0 ) {

        // Obtain key and value
        // Commands are delimited with a colon
        size_t colon_pos = line.find(":");
        std::string key = line.substr(0, colon_pos);
        key = key.substr(key.find_first_not_of(" "), key.find_last_not_of(" ")+1);
        std::string value = line.substr(colon_pos+1, line.size());
        value = value.substr(value.find_first_not_of(" "), value.find_last_not_of(" ")+1);

        if (key == "genesis input file") {
          genesis_file_name_ = value;
        }
        else if (key == "rve genesis input file") {
          rve_genesis_file_name_ = value;
        }
        else if (key == "exodus output file") {
          exodus_file_name_ = value;
        }
        else if (key == "use two level mesh decomposition") {
          std::string value_upper_case(value);
          std::transform(value_upper_case.begin(), value_upper_case.end(), value_upper_case.begin(), (int (*)(int))std::toupper);
          if (value_upper_case == "TRUE" || value_upper_case == "YES" || value_upper_case == "ON") {
            use_two_level_mesh_decomposition_ = true;
          }
          else if (value_upper_case == "FALSE" || value_upper_case == "NO" || value_upper_case == "OFF") {
            use_two_level_mesh_decomposition_ = false;
          }
          else {
            std::string msg = "\n**** Error in Parser::ReadFile(), unexpected value for \"use two level mesh decomposition\" " + value + "\n";
            throw std::logic_error(msg);
          }
        }
        else if (key == "write timing data file") {
          std::string value_upper_case(value);
          std::transform(value_upper_case.begin(), value_upper_case.end(), value_upper_case.begin(), (int (*)(int))std::toupper);
          if (value_upper_case == "TRUE" || value_upper_case == "YES" || value_upper_case == "ON") {
            write_timing_data_file_ = true;
          }
          else if (value_upper_case == "FALSE" || value_upper_case == "NO" || value_upper_case == "OFF") {
            write_timing_data_file_ = false;
          }
          else {
            std::string msg = "\n**** Error in Parser::ReadFile(), unexpected value for \"write timing data file\" " + value + "\n";
            throw std::logic_error(msg);
          }
        }
        else if (key == "time integration scheme") {
          time_integration_scheme_ = value;
        }
        else if (key == "nonlinear solver relative tolerance") {
          nonlinear_solver_relative_tolerance_ = std::atof(value.c_str());
        }
        else if (key == "nonlinear solver maximum iterations") {
          nonlinear_solver_max_iterations_ = std::atoi(value.c_str());
        }
        else if (key == "final time") {
          final_time_ = std::atof(value.c_str());
        }
        else if (key == "number of load steps") {
          num_load_steps_ = std::atoi(value.c_str());
        }
        else if (key == "output frequency") {
          output_frequency_ = std::atoi(value.c_str());
        }
        else if (key == "contact") {
          contact_string_ = value;
        }
        else if (key == "contact visualization") {

          std::stringstream ss(value);
          std::string val;
          std::vector<std::string> vals;
          while (ss >> val) {
            vals.push_back(val);
          }

          if (vals.size() != 6 ||
              vals[0] != "visualize_contact_entities" ||
              (vals[1] != "on" && vals[1] != "off") ||
              vals[2] != "visualize_bounding_boxes" ||
              (vals[3] != "on" && vals[3] != "off") ||
              vals[4] != "file_name") {
            std::string msg = "\n**** Error in Parser::ReadFile(), unexpected value for \"contact visualization\"\n";
            msg += "**** Allowable syntax is \"visualize_contatct_entities <on/off> visualize_bounding_boxes <on/off> file_name <file_name.e>\"\n";
            throw std::logic_error(msg);
          }

          visualize_contact_entities_ = false;
          if (vals[1] == "on") {
            visualize_contact_entities_ = true;
          }

          visualize_contact_bounding_boxes_ = false;
          if (vals[3] == "on") {
            visualize_contact_bounding_boxes_ = true;
          }

          contact_visualization_file_name_ = vals[5];
        }
        else if (key == "microscale output element ids") {
          std::stringstream ss(value);
          int global_id;
          while (ss >> global_id) {
            microscale_output_element_ids_.push_back(global_id);
          }
        }
        else if (key == "material parameters" ) {
          size_t space_pos = value.find(" ");
          std::string material_key = value.substr(0, space_pos);
          std::string material_props = value.substr(space_pos+1, value.size());
          material_strings_[material_key] = material_props;
        }
        else if (key == "macroscale block") {
          BlockProperties block_props(value);
          macroscale_blocks_[block_props.block_id_] = block_props;
        }
        else if (key == "microscale block") {
          BlockProperties block_props(value);
          microscale_blocks_[block_props.block_id_] = block_props;
        }
        else if (key == "microscale boundary condition strategy") {
          microscale_boundary_condition_strategy_ = value;
        }
        else if (key == "boundary condition") {
          boundary_condition_strings_.push_back(value);
        }
        else if (key == "output fields") {
          output_field_string_ = value;
        }
#ifdef NIMBLE_HAVE_UQ
	else if (key == "uq parameters") {
          size_t space_pos = value.find(" ");
          std::string material_key = value.substr(0, space_pos);
          std::string uq_params = value.substr(space_pos+1, value.size());
          uq_parameters_strings_[material_key] = uq_params;
        }
        else if (key == "uq model") {
          uq_model_string_ = value;
        }
#endif
        else{
          std::string msg = "\n**** Error in Parser::ReadFile(), unknown key " + key + "\n";
          throw std::logic_error(msg);
        }
      }
    }
  }
}
