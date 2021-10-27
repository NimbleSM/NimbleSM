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

#ifndef NIMBLE_MAIN_H
#define NIMBLE_MAIN_H

#include <memory>
#include <string>

#include "nimble_defs.h"
#include "nimble_cli.h"

namespace nimble {

class BlockMaterialInterfaceFactoryBase;
class ContactInterface;
class MaterialFactoryBase;
class Parser;

}  // namespace nimble

namespace nimble {

/// \brief Initialization routine that parses the input parameters
///
/// \param argc  Number of parameters on the command line
/// \param argv  List of parameters on the command line
/// \param parser  Parsing variable
///
/// \note The variable 'parser' should not be modified after calling
///  this function.
///
void
NimbleInitializeAndGetInput(int argc, char** argv, nimble::Parser& parser);

int
NimbleMain(
    const std::shared_ptr<nimble::MaterialFactoryBase>&               material_factory,
    std::shared_ptr<nimble::ContactInterface>                         contact_interface,
    const std::shared_ptr<nimble::BlockMaterialInterfaceFactoryBase>& block_material,
    const nimble::Parser&                                             parser);

void
NimbleFinalize(const nimble::EnvironmentFlags& env_flags);

}  // namespace nimble

#endif
