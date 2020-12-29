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

#include <memory>
#include <string>
#include <nimble_contact_interface.h>
#include "nimble_material_factory.h"

#include "nimble_mpi.h"

#include "nimble_parser.h"
#ifdef NIMBLE_HAVE_VT
#include <vt/transport.h>
#endif

int main(int argc, char *argv[]) {

  auto init_data = nimble::NimbleInitializeAndGetInput(argc, argv);

  int status = 0;
  {
    // Initialize VT if we need to
#ifdef NIMBLE_HAVE_VT
    MPI_Comm vt_comm = MPI_COMM_WORLD;
    auto vt_rt = ::vt::CollectiveOps::initialize(argc, argv, ::vt::no_workers, true, &vt_comm );
#endif
    std::shared_ptr<nimble::ContactInterface> contact_interface(new nimble::ContactInterface);
    std::shared_ptr<nimble::MaterialFactory> material_factory(new nimble::MaterialFactory);
    std::shared_ptr<nimble::Parser> parser(new nimble::Parser);

    status = nimble::NimbleMain(material_factory, contact_interface, parser,
                                init_data);
#ifdef NIMBLE_HAVE_VT
    while ( !vt_rt->isTerminated() )
      ::vt::runScheduler();
#endif
  }

  nimble::NimbleFinalize();

  return status;
}

