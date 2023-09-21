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

#ifndef NIMBLE_CLI_H
#define NIMBLE_CLI_H

#include <string>
#include <memory>

namespace nimble {
class CommandLineConfiguration
{
 public:

  CommandLineConfiguration(int argc, char** argv);
  CommandLineConfiguration(const CommandLineConfiguration &) = delete;
  CommandLineConfiguration(CommandLineConfiguration &&) noexcept;
  virtual ~CommandLineConfiguration();

  CommandLineConfiguration &operator=(const CommandLineConfiguration &) = delete;
  CommandLineConfiguration &operator=(CommandLineConfiguration &&) noexcept;

  virtual void ConfigureCommandLineArguments();

  int ParseAndGetErrorCode();

  bool UseTpetra() const noexcept;
  bool UseVT() const noexcept;
  bool UseUQ() const noexcept;

  const std::string &InputFilename() const noexcept;

  int &ArgC() noexcept;
  char **&ArgV() noexcept;

 private:

  struct impl;
  std::unique_ptr< impl > impl_;
};
}  // namespace nimble

#endif  // NIMBLE_CLI_H
