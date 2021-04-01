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

#include "nimble_contact_interface.h"
#include "nimble_kokkos.h"
#include "nimble_kokkos_block_material_interface_factory.h"
#include "nimble_kokkos_material_factory.h"

int
main(int argc, char** argv)
{
  nimble::Parser myParser;
  myParser.SetToUseKokkos();
  //
  // note: If the command line does not specify any parameter '--use-kokkos',
  // '--use-tpetra', or '--use-vt', we could choose to specify one of them with
  // either myParser.SetToUseKokkos() or myParser.SetToUseTpetra() or
  // myParser.SetToUseVT()
  //
  nimble::NimbleKokkosInitializeAndGetInput(argc, argv, myParser);

  {
    std::shared_ptr<MaterialFactoryType>                       material_factory(new nimble_kokkos::MaterialFactory);
    std::shared_ptr<nimble::BlockMaterialInterfaceFactoryBase> block_material_interface_factory(
        new nimble_kokkos::BlockMaterialInterfaceFactory);
    std::shared_ptr<nimble::ContactInterface> contact_interface = nullptr;
    if (myParser.HasContact())
      contact_interface = std::shared_ptr<nimble::ContactInterface>(new nimble::ContactInterface);

    nimble::NimbleKokkosMain(material_factory, contact_interface, block_material_interface_factory, myParser);
  }

  nimble::NimbleKokkosFinalize(myParser);

  return 0;
}
