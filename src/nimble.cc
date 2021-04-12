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

#include <exception>
#include <memory>

#include "nimble_block_material_interface_factory_base.h"
#include "nimble_contact_interface.h"
#include "nimble_main.h"
#include "nimble_material_factory.h"
#include "nimble_parser.h"

#ifdef NIMBLE_HAVE_KOKKOS
#include "nimble_kokkos_block_material_interface_factory.h"
#include "nimble_kokkos_material_factory.h"
#endif

int
main(int argc, char* argv[])
{
  //
  // Quick parsing of command line to check whether Kokkos, Tpetra, or VT is used.
  //
  nimble::EnvironmentFlags myFlags;
  nimble::parseCommandLine(argc, argv, myFlags);

  try {
    //
    // Create new Parser variable
    //
    nimble::Parser myParser;
    if (myFlags.use_kokkos_)
      myParser.SetToUseKokkos();
    else if (myFlags.use_tpetra_)
      myParser.SetToUseTpetra();
    else if (myFlags.use_vt_)
      myParser.SetToUseVT();

    nimble::NimbleInitializeAndGetInput(argc, argv, myParser);

    std::shared_ptr<nimble::MaterialFactoryBase> material_factory = nullptr;
#ifdef NIMBLE_HAVE_KOKKOS
    if (myParser.UseKokkos()) {
      material_factory = std::shared_ptr<nimble::MaterialFactoryBase>(new nimble_kokkos::MaterialFactory);
    } else
#endif
    {
      material_factory = std::shared_ptr<nimble::MaterialFactoryBase>(new nimble::MaterialFactory);
    }

    std::shared_ptr<nimble::BlockMaterialInterfaceFactoryBase> block_material = nullptr;
#ifdef NIMBLE_HAVE_KOKKOS
    if (myParser.UseKokkos()) {
      block_material =
          std::shared_ptr<nimble::BlockMaterialInterfaceFactoryBase>(new nimble_kokkos::BlockMaterialInterfaceFactory);
    }
#endif

    std::shared_ptr<nimble::ContactInterface> contact_interface = nullptr;
    if (myParser.HasContact())
      contact_interface = std::shared_ptr<nimble::ContactInterface>(new nimble::ContactInterface);

    int status = nimble::NimbleMain(material_factory, contact_interface, block_material, myParser);

  } 
  catch (std::exception &e) {
    std::cerr << "Standard exception: " << e.what() << std::endl;
  } 
  catch (...) {
    std::cerr << "\n !!! NimbleSM Simulation did not end correctly !!!\n\n";
  }

  nimble::NimbleFinalize(myFlags);

  return 0;
}
