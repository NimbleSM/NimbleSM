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

#ifndef NIMBLE_KOKKOS_H
#define NIMBLE_KOKKOS_H

#include <string>
#include <memory>

namespace nimble_kokkos { class MaterialFactory; }
namespace nimble { class ContactInterface; }
namespace nimble_kokkos { class BlockMaterialInterfaceFactory; }
namespace nimble { class Parser; }

namespace nimble {

struct NimbleKokkosInitData {
  std::string input_deck_name = "";
  int my_mpi_rank = -1;
  int num_mpi_ranks = -1;
};


NimbleKokkosInitData NimbleKokkosInitializeAndGetInput(int argc, char* argv[]);

void NimbleKokkosMain(std::shared_ptr<nimble_kokkos::MaterialFactory> material_factory,
                     std::shared_ptr<nimble::ContactInterface> contact_interface,
                     std::shared_ptr<nimble_kokkos::BlockMaterialInterfaceFactory> block_material_interface_factory,
                     std::shared_ptr<nimble::Parser> parser,
                     const NimbleKokkosInitData& init_data);

int NimbleKokkosFinalize(const NimbleKokkosInitData& init_data);

}

#endif
