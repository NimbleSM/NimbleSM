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

#include "nimble_timing_utils.h"

#include <fstream>
#include <iostream>
#include <sstream>

namespace nimble {

TimingInfo::TimingInfo(int i, std::string basicString,
                       double t_sim, double t_internal, double t_contact,
                       double t_exodus, double t_reduce)
    : num_ranks(i), time_stamp(std::move(basicString)),
      total_simulation_time(t_sim),
      total_internal_force_time(t_internal),
      total_contact_time(t_contact),
      total_exodus_write_time(t_exodus),
      total_vector_reduction_time(t_reduce)
    {}

void TimingInfo::BinaryWrite() const {

  std::stringstream timinginfo;
  timinginfo << num_ranks << "\t" << total_simulation_time
             << "\t" << total_internal_force_time
             << "\t" << total_contact_time
             << "\t" << total_exodus_write_time
             << "\t" << total_vector_reduction_time
             << "\n";
  std::string timingdatastr = timinginfo.str();

  std::stringstream filename_stream;
  filename_stream << "nimble_timing_data_n" << num_ranks << "_"
                  << time_stamp
                  << ".log";
  std::string timingfilenamestr = filename_stream.str();
  for (char &c : timingfilenamestr) {
    if (c == ' ')
      c = '-';
  }

  std::fstream fs(timingfilenamestr, fs.binary | fs.trunc | fs.in | fs.out);
  if (!fs.is_open()) {
    std::cerr << "Failed to open timing data file" 
              << " (" << timingfilenamestr << ") "
              << std::endl;
  } else {
    fs << timingdatastr;
    fs.flush();
    fs.close();
  }

}

} // namespace nimble

