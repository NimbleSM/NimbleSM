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

#include "nimble.h"
#include <exception>
#include <memory>

#ifdef NIMBLE_HAVE_MPI
#include <mpi.h>
#endif

#include "nimble_block_material_interface_factory_base.h"
#include "nimble_contact_interface.h"
#include "nimble_material_factory.h"
#include "nimble_parser.h"
#include "nimble_version.h"
#include "nimble_cli.h"
#include "nimble_genesis_mesh.h"
#include "nimble_data_manager.h"
#include "integrators/integrator_base.h"
#include "integrators/explicit_time_integrator.h"
#include "integrators/quasistatic_time_integrator.h"

#ifdef NIMBLE_HAVE_KOKKOS
#include "nimble_kokkos_block_material_interface_factory.h"
#include "nimble_kokkos_material_factory.h"
#endif

#ifdef NIMBLE_HAVE_VT
#include <vt/transport.h>
#endif

namespace nimble {
struct NimbleApplication::impl
{
#ifdef NIMBLE_HAVE_VT
  ::vt::RuntimePtrType vt_rt_ = nullptr;
#endif

  std::unique_ptr< Parser > parser_;
  std::unique_ptr< CommandLineConfiguration > cli_config_;
  int rank_ = 0;
  int num_ranks_ = 1;

  // Customizable interfaces and factories
  std::shared_ptr< MaterialFactoryBase > material_factory_;
  std::shared_ptr< BlockMaterialInterfaceFactoryBase > block_material_interface_factory_;
  std::shared_ptr< ContactInterface > contact_interface_;
  std::unique_ptr< IntegratorBase > integrator_;
};

NimbleApplication::NimbleApplication()
  : impl_{ std::make_unique< impl >() }
{

}

NimbleApplication::NimbleApplication( NimbleApplication && ) noexcept = default;

NimbleApplication::~NimbleApplication() = default;

NimbleApplication &
NimbleApplication::operator=( NimbleApplication && ) noexcept = default;

int
NimbleApplication::Run(int argc, char **argv) {
  impl_->cli_config_ = CreateCLI(argc, argv);

  impl_->cli_config_->ConfigureCommandLineArguments();
  int err = impl_->cli_config_->ParseAndGetErrorCode();
  if (err != 0)
    return err;

  impl_->parser_ = CreateParser();

  if (impl_->cli_config_->UseKokkos())
    impl_->parser_->SetToUseKokkos();
  if (impl_->cli_config_->UseTpetra())
    impl_->parser_->SetToUseTpetra();
  if (impl_->cli_config_->UseVT())
    impl_->parser_->SetToUseVT();
  if (impl_->cli_config_->UseUQ())
    impl_->parser_->SetToUseUQ();

  impl_->parser_->SetInputFilename(impl_->cli_config_->InputFilename());

  InitializeSubsystems();

  // Run main; if using VT this will run it in an epoch
  int status = 0;
  try {
    status = ExecMain();
  } catch (std::exception& e) {
    std::cerr << "Standard exception: " << e.what() << std::endl;
    return 1;
  } catch (...) {
    std::cerr << "\n !!! NimbleSM Simulation did not end correctly !!!\n\n";
    return 1;
  }

  // We're done, finalize
  FinalizeSubsystems();

  // Run is technically re-entrant so make sure contact interface is de-initialized
  // because if the aprser changes it could end up with strange results
  // We don't need to do it for the factories since they will be overwritten if
  // run() is called again
  impl_->contact_interface_ = {};

  return status;
}

int
NimbleApplication::Rank() const noexcept
{
  return impl_->rank_;
}

int
NimbleApplication::NumRanks() const noexcept
{
  return impl_->num_ranks_;
}

std::shared_ptr< ContactInterface >
NimbleApplication::GetContactInterface()
{
  return impl_->contact_interface_;
}

Parser &
NimbleApplication::GetParser()
{
  return *impl_->parser_;
}

void
NimbleApplication::InitializeSubsystems()
{
#ifdef NIMBLE_HAVE_TRILINOS
  if (parser.UseTpetra()) {
    auto sguard = new Tpetra::ScopeGuard(&argc, &argv);
    parser.ResetTpetraScope(sguard);
    auto comm = Tpetra::getDefaultComm();
    num_ranks_ = comm->getSize();
    rank_   = comm->getRank();
  } else
#endif
#ifdef NIMBLE_HAVE_MPI
{
    MPI_Init(&impl_->cli_config_->ArgC(), &impl_->cli_config_->ArgV());
    int mpi_err = MPI_Comm_rank(MPI_COMM_WORLD, &impl_->rank_);
    if (mpi_err != MPI_SUCCESS) { NIMBLE_ABORT("\nError:  MPI_Comm_rank() returned nonzero error code.\n"); }
    mpi_err = MPI_Comm_size(MPI_COMM_WORLD, &impl_->num_ranks_);
    if (mpi_err != MPI_SUCCESS) { NIMBLE_ABORT("\nError:  MPI_Comm_size() returned nonzero error code.\n"); }
  }
#endif

#ifdef NIMBLE_HAVE_KOKKOS
Kokkos::initialize(impl_->cli_config_->ArgC(), impl_->cli_config_->ArgV());
#endif

  // Banner
  if (impl_->rank_ == 0) {
    std::cout << "\n-- NimbleSM" << std::endl;
    std::cout << "-- version " << nimble::NimbleVersion() << "\n";
    if (impl_->parser_->UseKokkos()) {
      std::cout << "-- Using Kokkos interface \n";
    }
    if (impl_->parser_->UseTpetra()) {
      std::cout << "-- Using Tpetra interface \n";
    }
    if (impl_->parser_->UseVT()) {
      std::cout << "-- Using VT runtime \n";
    }
    std::cout << "-- Number of ranks = " << impl_->num_ranks_ << "\n";
    std::cout << std::endl;
  }

  // Initialize VT if needed
#ifdef NIMBLE_HAVE_VT
  if (impl_->parser_->UseVT()) {
    MPI_Comm vt_comm = MPI_COMM_WORLD;
    int argc = impl_->cli_config_->ArgC();
    char **argv = impl_->cli_config_->ArgV();
    impl_->vt_rt_    = ::vt::CollectiveOps::initialize(argc, argv, ::vt::no_workers, true, &vt_comm);
  }
#endif

  impl_->parser_->SetRankID(impl_->rank_);
  impl_->parser_->SetNumRanks(impl_->num_ranks_);

  impl_->parser_->Initialize();

  // Now initiate our factories and interfaces
  impl_->material_factory_ = CreateMaterialFactory();
  impl_->block_material_interface_factory_ = CreateBlockMaterialInterfaceFactory();
  if ( impl_->parser_->HasContact() )
    impl_->contact_interface_ = CreateContactInterface();
}

void NimbleApplication::FinalizeSubsystems()
{
#ifdef NIMBLE_HAVE_VT
  if (impl_->parser_->UseVT())
    ::vt::finalize(std::move(impl_->vt_rt_));
#endif

#ifdef NIMBLE_HAVE_KOKKOS
  if (impl_->parser_->UseKokkos())
    Kokkos::finalize();
#endif

#ifdef NIMBLE_HAVE_TRILINOS
  if (!impl_->parser_->UseTpetra) {
#ifdef NIMBLE_HAVE_MPI
    MPI_Finalize();
#endif
  }
#else
#ifdef NIMBLE_HAVE_MPI
  MPI_Finalize();
#endif
#endif
}

std::unique_ptr< Parser >
NimbleApplication::CreateParser()
{
  return std::make_unique< Parser >();
}

std::unique_ptr< CommandLineConfiguration >
NimbleApplication::CreateCLI(int argc, char **argv)
{
  return std::make_unique< CommandLineConfiguration >(argc, argv);
}

std::unique_ptr< MaterialFactoryBase >
NimbleApplication::CreateMaterialFactory()
{
#ifdef NIMBLE_HAVE_KOKKOS
  if (impl_->parser_->UseKokkos()) {
    return std::make_unique<nimble_kokkos::MaterialFactory>();
  } else
#endif
  {
    return std::make_unique<nimble::MaterialFactory>();
  }
}

std::unique_ptr< BlockMaterialInterfaceFactoryBase >
NimbleApplication::CreateBlockMaterialInterfaceFactory()
{
#ifdef NIMBLE_HAVE_KOKKOS
  if (impl_->parser_->UseKokkos()) {
    return std::make_unique<nimble_kokkos::BlockMaterialInterfaceFactory>();
  } else
#endif
  {
    return {};
  }
}

std::unique_ptr< ContactInterface >
NimbleApplication::CreateContactInterface()
{
  return std::make_unique< ContactInterface >();
}

std::unique_ptr< IntegratorBase >
NimbleApplication::CreateIntegrator(DataManager &data_manager,
                                    GenesisMesh &mesh)
{
  auto time_integration_scheme = impl_->parser_->TimeIntegrationScheme();
  if (time_integration_scheme == "explicit") {
    return std::make_unique< ExplicitTimeIntegrator >(*this, mesh, data_manager);
  } else if (time_integration_scheme == "quasistatic") {
    return std::make_unique< QuasistaticTimeIntegrator >(*this, mesh, data_manager);
  }

  return {};
}

int
NimbleApplication::ExecMain()
{
  int ret = 0;
#ifdef NIMBLE_HAVE_VT
  if (impl_->parser_->UseVT()) {
    vt::runInEpochCollective([this, &ret]() {
      ret = MainLoop();
    });
  } else
#endif
  {
    ret = MainLoop();
  }

  return ret;
}

int
NimbleApplication::MainLoop()
{
  const int my_rank   = impl_->rank_;
  const int num_ranks = impl_->num_ranks_;

  // Read the mesh
  nimble::GenesisMesh mesh;
  {
    std::string genesis_file_name = nimble::IOFileName(impl_->parser_->GenesisFileName(), "g", "", my_rank, num_ranks);
    mesh.ReadFile(genesis_file_name);
  }

  std::string tag                = "out";
  std::string output_exodus_name = nimble::IOFileName(impl_->parser_->ExodusFileName(), "e", tag, my_rank, num_ranks);

  int dim       = mesh.GetDim();
  int num_nodes = static_cast<int>(mesh.GetNumNodes());

  nimble::DataManager data_manager(*impl_->parser_, mesh);
  data_manager.SetBlockMaterialInterfaceFactory(impl_->block_material_interface_factory_);

  if (my_rank == 0) {
    std::cout << "\n";
    if (num_ranks == 1) {
      std::cout << " Number of Nodes = " << num_nodes << "\n";
      std::cout << " Number of Elements = " << mesh.GetNumElements() << "\n";
    }
    std::cout << " Number of Global Blocks = " << mesh.GetNumGlobalBlocks() << "\n";
    std::cout << "\n";
    std::cout << " Number of Ranks         = " << num_ranks << "\n";
#ifdef _OPENMP
    std::cout << " Number of Threads       = " << omp_get_max_threads() << "\n";
#endif
    std::cout << "\n";
  }

  auto model_data = data_manager.GetModelData();
  model_data->InitializeBlocks(data_manager, impl_->material_factory_);

  //
  // Initialize the output file
  //

  data_manager.InitializeOutput(output_exodus_name);

  auto integrator = CreateIntegrator(data_manager, mesh);

  return integrator->Integrate();
}
}
