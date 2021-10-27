/*
// @HEADER
//  ************************************************************************
//
//                                 NimbleSM
//                              Copyright 2018
//    National Technology & Engineering Solutions of Sandia, LLC (NTESS)
//
//  Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government
//  retains certain rights in this software.
//
//  Redistribution and use in source and binary forms, with or without
//  modification, are permitted provided that the following conditions are
//  met:
//
//  1. Redistributions of source code must retain the above copyright
//  notice, this list of conditions and the following disclaimer.
//
//  2. Redistributions in binary form must reproduce the above copyright
//  notice, this list of conditions and the following disclaimer in the
//  documentation and/or other materials provided with the distribution.
//
//  3. Neither the name of the Corporation nor the names of the
//  contributors may be used to endorse or promote products derived from
//  this software without specific prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY EXPRESS OR IMPLIED
//  WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
//  MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN
//  NO EVENT SHALL NTESS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
//  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
//  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
//  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
//  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
//  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
//  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
//  Questions?  Contact David Littlewood (djlittl@sandia.gov)
//
//  ************************************************************************
// @HEADER
*/

#include "nimble_cli.h"

#include <CLI11/CLI11.hpp>

namespace nimble {
int
parseCommandLine(int argc, char** argv, nimble::EnvironmentFlags& myFlags)
{
  CLI::App nimble_app{"A Lagrangian finite-element code for solid mechanics."};
  nimble_app.allow_extras();

  auto* kokkos_opt = nimble_app.add_flag("--use_kokkos", myFlags.use_kokkos_, "Use Kokkos for parallel compute")
                         ->check([](const std::string&) -> std::string {
#ifdef NIMBLE_HAVE_KOKKOS
                           return {};  // Empty string on success
#else
                           return "NimbleSM was not built with Kokkos support";
#endif
                         });

  auto* tpetra_opt = nimble_app.add_flag("--use_tpetra", myFlags.use_tpetra_, "Use TPetra solvers")
                         ->check([](const std::string&) -> std::string {
#ifdef NIMBLE_HAVE_TRILINOS
                           return {};  // Empty string on success
#else
                           return "NimbleSM was not built with Trilinos support";
#endif
                         });

  nimble_app.add_flag("--use_vt", myFlags.use_vt_, "Use DARMA/vt for asynchronous distributed contact tasking")
          ->check([](const std::string&) -> std::string {
#ifdef NIMBLE_HAVE_VT
            return {};  // Empty string on success
#else
            return "NimbleSM was not built with DARMA/vt support";
#endif
          });

  kokkos_opt->excludes(tpetra_opt);
  tpetra_opt->excludes(kokkos_opt);

  try {
    nimble_app.parse(argc, argv);
  } catch (const CLI::ParseError& e) {
    return nimble_app.exit(e);
  }
  return 0;
}
}  // namespace nimble