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

#ifndef NIMBLE_OUTPUT_EXODUS_H
#define NIMBLE_OUTPUT_EXODUS_H

#include "nimble_genesis_mesh.h"

#ifdef NIMBLE_HAVE_DARMA
  #include "darma.h"
#else
  #include <vector>
  #include <map>
#endif

namespace nimble {

  class ExodusOutput {

  public:

    ExodusOutput() : filename_("none"), CPU_word_size_(sizeof(double)), IO_word_size_(sizeof(double)),
      dim_(0), num_nodes_(0), num_elements_(0), num_blocks_(0),
      num_node_sets_(0), num_side_sets_(0), exodus_write_count_(0) {}

#ifdef NIMBLE_HAVE_DARMA
    template<typename ArchiveType>
    void serialize(ArchiveType& ar) {
			ar | filename_ | CPU_word_size_ | IO_word_size_ | dim_ | num_nodes_ | num_elements_ | num_blocks_ | block_ids_ | num_node_sets_ | num_side_sets_ | exodus_write_count_ | elem_data_index_;
		}
#endif

    void Initialize(std::string const & filename,
                    GenesisMesh const & genesis_mesh);

    virtual ~ExodusOutput() {}

    std::string GetFileName() const { return filename_; }

    void InitializeDatabase(GenesisMesh const & genesis_mesh,
			    std::vector<std::string> const & global_data_names,
			    std::vector<std::string> const & node_data_names,
          std::map< int, std::vector<std::string> > const & elem_data_names,
          std::map< int, std::vector<std::string> > const & derived_elem_data_names);

    void WriteStep(double time,
		   std::vector<double> const & global_data,
		   std::vector< std::vector<double> > const & node_data,
       std::map< int, std::vector< std::string > > const & elem_data_names,
       std::map< int, std::vector< std::vector<double> > > const & elem_data,
       std::map< int, std::vector< std::string > > const & derived_elem_data_names,
       std::map< int, std::vector< std::vector<double> > > const & derived_elem_data);

    void InitializeDatabaseTextFile(GenesisMesh const & genesis_mesh,
       std::vector<std::string> const & global_data_names,
       std::vector<std::string> const & node_data_names,
       std::map< int, std::vector<std::string> > const & elem_data_names,
       std::map< int, std::vector<std::string> > const & derived_elem_data_names);

    void WriteStepTextFile(double time,
		   std::vector<double> const & global_data,
		   std::vector< std::vector<double> > const & node_data,
       std::map< int, std::vector< std::string > > const & elem_data_names,
       std::map< int, std::vector< std::vector<double> > > const & elem_data,
       std::map< int, std::vector< std::string > > const & derived_elem_data_names,
       std::map< int, std::vector< std::vector<double> > > const & derived_elem_data);

  protected:

    void ReportExodusError(int error_code,
			   const char *method_name,
			   const char *exodus_method_name);

    void WriteQARecord(int exodus_file_id);

    std::string filename_;
    int CPU_word_size_;
    int IO_word_size_;
    int dim_;
    int num_nodes_;
    int num_elements_;
    int num_blocks_;
    std::vector<int> block_ids_;
    int num_node_sets_;
    int num_side_sets_;
    int exodus_write_count_;
    std::map<std::string, int> elem_data_index_;
  };

} // namespace nimble

#endif // NIMBLE_OUTPUT_EXODUS_H
