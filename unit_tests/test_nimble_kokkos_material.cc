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

#include <gtest/gtest.h>
#include <nimble_kokkos_defs.h>
#include <nimble_kokkos_material_factory.h>
#include <nimble_material.h>

#include <iostream>
#include <string>

using namespace nimble_kokkos;

inline void
execute_material_on_device(nimble::Material* mat)
{
  DeviceFullTensorIntPtView all_dn("d", 1);
  DeviceFullTensorIntPtView all_dnp1("dnp1", 1);
  DeviceSymTensorIntPtView  all_sn("sn", 1);
  DeviceSymTensorIntPtView  all_snp1("snp1", 1);
  auto                      dn   = Kokkos::subview(all_dn, 0, 0, Kokkos::ALL());
  auto                      dnp1 = Kokkos::subview(all_dnp1, 0, 0, Kokkos::ALL());
  auto                      sn   = Kokkos::subview(all_sn, 0, 0, Kokkos::ALL());
  auto                      snp1 = Kokkos::subview(all_snp1, 0, 0, Kokkos::ALL());

  Kokkos::parallel_for(
      1, KOKKOS_LAMBDA(const int i) {
        auto dens = mat->GetDensity();
        printf("Density = %f\n", dens);
        mat->GetStress(0., 0., dn, dnp1, sn, snp1);
      });
}

TEST(nimble_kokkos_material, can_create_on_host_and_device)
{
  const std::string material_string = "neohookean bulk_modulus 1.0e6 shear_modulus 5.e5 density 1.0e3";
  MaterialFactory   kokkosFactory;
  kokkosFactory.parse_and_create(material_string, 0);

  auto mat_dev = kokkosFactory.get_material_device();
  ASSERT_NE(nullptr, mat_dev);

  ASSERT_NO_THROW(execute_material_on_device(mat_dev));

  Kokkos::kokkos_free(mat_dev);
}

class TestKokkosMaterialFactory : public MaterialFactory
{
 public:
  TestKokkosMaterialFactory() : MaterialFactory() {}
  ~TestKokkosMaterialFactory() = default;

  std::shared_ptr<nimble::MaterialParameters>
  parse_string(const std::string& params) const
  {
    return ParseMaterialParametersString(params);
  }

  void
  create() override
  {
  }
};

template <typename MatType>
struct TestMatAlloc
{
  TestMatAlloc(MatType& host_mat)
      : mat_device(static_cast<MatType*>(Kokkos::kokkos_malloc<>("Material", sizeof(MatType))))
  {
    allocate(host_mat, mat_device);
  }

  void
  allocate(MatType& host_mat, MatType* dev_mat)
  {
    assert(dev_mat != nullptr);
    Kokkos::parallel_for(
        1, KOKKOS_LAMBDA(int) { new (dev_mat) MatType(host_mat); });
    Kokkos::fence();
  }

  ~TestMatAlloc() { Kokkos::kokkos_free(mat_device); }

  MatType* mat_device;
};

TEST(nimble_kokkos_material, can_create_and_run_directly)
{
  using MatType                            = nimble::ElasticMaterial;
  const std::string        material_string = "neohookean bulk_modulus 1.0e6 shear_modulus 5.e5 density 1.0e3";
  auto                     params          = TestKokkosMaterialFactory().parse_string(material_string);
  std::shared_ptr<MatType> mat_host(new MatType(*params));

  TestMatAlloc<MatType> alloc(*mat_host);
  nimble::Material*     mat_abstract = alloc.mat_device;
  ASSERT_NO_THROW(execute_material_on_device(mat_abstract));
}
