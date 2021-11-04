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

#include "nimble_parser.h"

#ifdef NIMBLE_HAVE_VT
// From VT, if not we get some nasty conflicting symbols
#include <CLI/CLI11.hpp>
#else
#include <nimble_tpl/CLI11/CLI11.hpp>
#endif

namespace nimble {

struct CommandLineConfiguration::impl
{
  CLI::App cli_app_{"A Lagrangian finite-element code for solid mechanics."};
  int argc_ = 0;
  char **argv_ = nullptr;

  bool use_kokkos_ = false;
  bool use_tpetra_ = false;
  bool use_vt_ = false;
  bool use_uq_ = false;
  std::string input_filename_;
};

CommandLineConfiguration::CommandLineConfiguration(int argc, char** argv)
  : impl_{std::make_unique< impl >()}
{
  impl_->argc_ = argc;
  impl_->argv_ = argv;
}

CommandLineConfiguration::CommandLineConfiguration(CommandLineConfiguration &&) noexcept = default;

CommandLineConfiguration::~CommandLineConfiguration() = default;

CommandLineConfiguration &
CommandLineConfiguration::operator=(CommandLineConfiguration &&) noexcept = default;

void
CommandLineConfiguration::ConfigureCommandLineArguments() {
  impl_->cli_app_.allow_extras();

  auto* kokkos_opt = impl_->cli_app_.add_flag("--use_kokkos", impl_->use_kokkos_, "Use Kokkos for parallel compute")
                         ->check([](const std::string&) -> std::string {
#ifdef NIMBLE_HAVE_KOKKOS
                           return {};  // Empty string on success
#else
                           return "NimbleSM was not built with Kokkos support";
#endif
                         });

  auto* tpetra_opt = impl_->cli_app_.add_flag("--use_tpetra", impl_->use_tpetra_, "Use TPetra solvers")
                         ->check([](const std::string&) -> std::string {
#ifdef NIMBLE_HAVE_TRILINOS
                           return {};  // Empty string on success
#else
                           return "NimbleSM was not built with Trilinos support";
#endif
                         });

  impl_->cli_app_.add_flag("--use_vt", impl_->use_vt_, "Use DARMA/vt for asynchronous distributed contact tasking")
          ->check([](const std::string&) -> std::string {
#ifdef NIMBLE_HAVE_VT
            return {};  // Empty string on success
#else
            return "NimbleSM was not built with DARMA/vt support";
#endif
          });

  impl_->cli_app_.add_flag("--use_uq", impl_->use_vt_, "Use UQ capability")->check([](const std::string&) -> std::string {
#ifdef NIMBLE_HAVE_UQ
    return {};  // Empty string on success
#else
    return "NimbleSM was not built with UQ support";
#endif
  });

  impl_->cli_app_.add_option( "INPUT_FILE", impl_->input_filename_)->required()->check(CLI::ExistingFile);

  kokkos_opt->excludes(tpetra_opt);
  tpetra_opt->excludes(kokkos_opt);
}

int
CommandLineConfiguration::ParseAndGetErrorCode() {
  try {
    impl_->cli_app_.parse(impl_->argc_, impl_->argv_);
  } catch (const CLI::ParseError& e) {
    return impl_->cli_app_.exit(e);
  }

  return 0;
}

bool
CommandLineConfiguration::UseKokkos() const noexcept
{
  return impl_->use_kokkos_;
}

bool
CommandLineConfiguration::UseTpetra() const noexcept
{
  return impl_->use_tpetra_;
}

bool
CommandLineConfiguration::UseVT() const noexcept
{
  return impl_->use_vt_;
}

bool
CommandLineConfiguration::UseUQ() const noexcept
{
  return impl_->use_uq_;
}

const std::string &
CommandLineConfiguration::InputFilename() const noexcept
{
  return impl_->input_filename_;
}

int &
CommandLineConfiguration::ArgC() noexcept
{
  return impl_->argc_;
}

char **&
CommandLineConfiguration::ArgV() noexcept
{
  return impl_->argv_;
}

}  // namespace nimble