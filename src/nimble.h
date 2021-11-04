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

#ifndef NIMBLE_H
#define NIMBLE_H

#include <memory>

namespace nimble {
class Parser;
class CommandLineConfiguration;
class MaterialFactoryBase;
class BlockMaterialInterfaceFactoryBase;
class ContactInterface;
class IntegratorBase;
class DataManager;
class GenesisMesh;

class NimbleApplication
{
 public:
  NimbleApplication();
  NimbleApplication( const NimbleApplication & ) = delete;
  NimbleApplication( NimbleApplication && ) noexcept;
  virtual ~NimbleApplication();

  NimbleApplication &operator=( const NimbleApplication & ) = delete;
  NimbleApplication &operator=( NimbleApplication && ) noexcept;

  int Run(int argc, char **argv);

  int Rank() const noexcept;
  int NumRanks() const noexcept;

  std::shared_ptr< ContactInterface > GetContactInterface();
  Parser &GetParser();

 private:

  struct impl;

  virtual std::unique_ptr< Parser > CreateParser();
  virtual std::unique_ptr< CommandLineConfiguration > CreateCLI(int argc, char **argv);
  virtual std::unique_ptr< MaterialFactoryBase > CreateMaterialFactory();
  virtual std::unique_ptr< BlockMaterialInterfaceFactoryBase > CreateBlockMaterialInterfaceFactory();
  virtual std::unique_ptr< ContactInterface > CreateContactInterface();
  virtual std::unique_ptr< IntegratorBase > CreateIntegrator(DataManager &data_manager,
                                                             GenesisMesh &mesh);

  void InitializeSubsystems();
  void FinalizeSubsystems();
  int ExecMain();
  int MainLoop();

  std::unique_ptr< impl > impl_;
};
}  // namespace nimble

#endif  // NIMBLE_H
