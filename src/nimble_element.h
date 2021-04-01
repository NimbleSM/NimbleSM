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

#ifndef NIMBLE_ELEMENT_H
#define NIMBLE_ELEMENT_H

#include "nimble_defs.h"

namespace nimble {

class Element
{
 public:
  NIMBLE_FUNCTION
  Element() {}

  NIMBLE_FUNCTION
  virtual ~Element() {}

  virtual int
  Dim() = 0;

  virtual int
  NumNodesPerElement() = 0;

  virtual int
  NumIntegrationPointsPerElement() = 0;

  virtual void
  ComputeLumpedMass(const double density, const double* node_reference_coords, double* lumped_mass) const = 0;

#ifdef NIMBLE_HAVE_KOKKOS
  NIMBLE_FUNCTION
  virtual void
  ComputeLumpedMass(
      const double                                   density,
      nimble_kokkos::DeviceVectorNodeGatheredSubView node_reference_coords,
      nimble_kokkos::DeviceScalarNodeGatheredSubView lumped_mass) const = 0;
#endif

  virtual double
  ComputeCharacteristicLength(const double* node_coords) = 0;

  virtual void
  ComputeVolumeAverage(
      const double* node_current_coords,
      int           num_quantities,
      const double* int_pt_quantities,
      double&       volume,
      double*       volume_averaged_quantity) const = 0;

#ifdef NIMBLE_HAVE_KOKKOS

  NIMBLE_FUNCTION
  virtual void
  ComputeVolume(
      nimble_kokkos::DeviceVectorNodeGatheredSubView node_reference_coords,
      nimble_kokkos::DeviceVectorNodeGatheredSubView node_displacements,
      nimble_kokkos::DeviceScalarElemSingleEntryView elem_volume) const = 0;

  NIMBLE_FUNCTION
  virtual void
  ComputeVolumeAverageSymTensor(
      nimble_kokkos::DeviceVectorNodeGatheredSubView    node_reference_coords,
      nimble_kokkos::DeviceVectorNodeGatheredSubView    node_displacements,
      nimble_kokkos::DeviceSymTensorIntPtSubView        int_pt_quantities,
      nimble_kokkos::DeviceSymTensorElemSingleEntryView vol_ave_quantity) const = 0;

  NIMBLE_FUNCTION
  virtual void
  ComputeVolumeAverageFullTensor(
      nimble_kokkos::DeviceVectorNodeGatheredSubView     node_reference_coords,
      nimble_kokkos::DeviceVectorNodeGatheredSubView     node_displacements,
      nimble_kokkos::DeviceFullTensorIntPtSubView        int_pt_quantities,
      nimble_kokkos::DeviceFullTensorElemSingleEntryView vol_ave_quantity) const = 0;
#endif

  virtual void
  ComputeDeformationGradients(
      const double* node_reference_coords,
      const double* node_current_coords,
      double*       deformation_gradients) const = 0;

#ifdef NIMBLE_HAVE_KOKKOS
  NIMBLE_FUNCTION
  virtual void
  ComputeDeformationGradients(
      nimble_kokkos::DeviceVectorNodeGatheredSubView node_reference_coords,
      nimble_kokkos::DeviceVectorNodeGatheredSubView node_displacements,
      nimble_kokkos::DeviceFullTensorIntPtSubView    deformation_gradients) const = 0;
#endif

  virtual void
  ComputeTangent(const double* node_reference_coords, const double* node_current_coords, double* tangent) = 0;

  virtual void
  ComputeNodalForces(const double* node_current_coords, const double* int_pt_stresses, double* node_forces) = 0;

#ifdef NIMBLE_HAVE_KOKKOS
  NIMBLE_FUNCTION
  virtual void
  ComputeNodalForces(
      nimble_kokkos::DeviceVectorNodeGatheredSubView node_reference_coords,
      nimble_kokkos::DeviceVectorNodeGatheredSubView node_displacements,
      nimble_kokkos::DeviceSymTensorIntPtSubView     element_stress_step_np1_d,
      nimble_kokkos::DeviceVectorNodeGatheredSubView element_internal_force_d) const = 0;
#endif

  NIMBLE_INLINE_FUNCTION
  double
  Invert3x3(double mat[][3], double inv[][3]) const;

  NIMBLE_INLINE_FUNCTION
  void
  LU_Decompose(double mat[][3], int index[]) const;

  NIMBLE_INLINE_FUNCTION
  void
  LU_Solve(double a[][3], int index[3], double b[3]) const;

  NIMBLE_INLINE_FUNCTION
  void
  LU_Invert(double mat[][3], double inv[][3]) const;

  NIMBLE_INLINE_FUNCTION
  double
  MatrixInverseCheckCorrectness(double mat[][3], double inv[][3]) const;

 protected:
  const int K_S_XX_ = 0;
  const int K_S_YY_ = 1;
  const int K_S_ZZ_ = 2;
  const int K_S_XY_ = 3;
  const int K_S_YZ_ = 4;
  const int K_S_ZX_ = 5;
  const int K_S_YX_ = 3;
  const int K_S_ZY_ = 4;
  const int K_S_XZ_ = 5;

  const int K_F_XX_ = 0;
  const int K_F_YY_ = 1;
  const int K_F_ZZ_ = 2;
  const int K_F_XY_ = 3;
  const int K_F_YZ_ = 4;
  const int K_F_ZX_ = 5;
  const int K_F_YX_ = 6;
  const int K_F_ZY_ = 7;
  const int K_F_XZ_ = 8;
};

class HexElement : public Element
{
 public:
  NIMBLE_FUNCTION
  HexElement();

  NIMBLE_FUNCTION
  virtual ~HexElement() {}

  int
  Dim()
  {
    return 3;
  }

  int
  NumNodesPerElement()
  {
    return 8;
  }

  int
  NumIntegrationPointsPerElement()
  {
    return 8;
  }

  void
  ComputeLumpedMass(const double density, const double* node_reference_coords, double* lumped_mass) const;

#ifdef NIMBLE_HAVE_KOKKOS
  NIMBLE_FUNCTION
  void
  ComputeLumpedMass(
      const double                                   density,
      nimble_kokkos::DeviceVectorNodeGatheredSubView node_reference_coords,
      nimble_kokkos::DeviceScalarNodeGatheredSubView lumped_mass) const;
#endif

  double
  ComputeCharacteristicLength(const double* node_coords);

  void
  ComputeVolumeAverage(
      const double* node_current_coords,
      int           num_quantities,
      const double* int_pt_quantities,
      double&       volume,
      double*       volume_averaged_quantity) const;

#ifdef NIMBLE_HAVE_KOKKOS

  NIMBLE_FUNCTION
  void
  ComputeVolume(
      nimble_kokkos::DeviceVectorNodeGatheredSubView node_reference_coords,
      nimble_kokkos::DeviceVectorNodeGatheredSubView node_displacements,
      nimble_kokkos::DeviceScalarElemSingleEntryView elem_volume) const;

  NIMBLE_FUNCTION
  void
  ComputeVolumeAverageSymTensor(
      nimble_kokkos::DeviceVectorNodeGatheredSubView    node_reference_coords,
      nimble_kokkos::DeviceVectorNodeGatheredSubView    node_displacements,
      nimble_kokkos::DeviceSymTensorIntPtSubView        int_pt_quantities,
      nimble_kokkos::DeviceSymTensorElemSingleEntryView vol_ave_quantity) const;

  NIMBLE_FUNCTION
  void
  ComputeVolumeAverageFullTensor(
      nimble_kokkos::DeviceVectorNodeGatheredSubView     node_reference_coords,
      nimble_kokkos::DeviceVectorNodeGatheredSubView     node_displacements,
      nimble_kokkos::DeviceFullTensorIntPtSubView        int_pt_quantities,
      nimble_kokkos::DeviceFullTensorElemSingleEntryView vol_ave_quantity) const;
#endif

  void
  ComputeDeformationGradients(
      const double* node_reference_coords,
      const double* node_current_coords,
      double*       deformation_gradients) const;

#ifdef NIMBLE_HAVE_KOKKOS
  NIMBLE_FUNCTION
  void
  ComputeDeformationGradients(
      nimble_kokkos::DeviceVectorNodeGatheredSubView node_reference_coords,
      nimble_kokkos::DeviceVectorNodeGatheredSubView node_displacements,
      nimble_kokkos::DeviceFullTensorIntPtSubView    deformation_gradients) const;
#endif

  void
  ComputeTangent(const double* node_current_coords, const double* material_tangent, double* element_tangent);

  void
  ComputeNodalForces(const double* node_current_coords, const double* int_pt_stresses, double* node_forces);

#ifdef NIMBLE_HAVE_KOKKOS
  NIMBLE_FUNCTION
  void
  ComputeNodalForces(
      nimble_kokkos::DeviceVectorNodeGatheredSubView node_reference_coords,
      nimble_kokkos::DeviceVectorNodeGatheredSubView node_displacements,
      nimble_kokkos::DeviceSymTensorIntPtSubView     int_pt_stresses,
      nimble_kokkos::DeviceVectorNodeGatheredSubView node_forces) const;
#endif

 protected:
  NIMBLE_FUNCTION
  void
  ShapeFunctionValues(const double* natural_coords, double* shape_function_values);

  NIMBLE_FUNCTION
  void
  ShapeFunctionDerivatives(const double* natural_coords, double* shape_function_derivatives);

 private:
  static constexpr int dim_         = 3;
  static constexpr int num_nodes_   = 8;
  static constexpr int num_int_pts_ = 8;

  double int_pts_[num_int_pts_ * dim_];
  double int_wts_[num_int_pts_];
  double shape_fcn_vals_[num_nodes_ * num_int_pts_];
  double shape_fcn_deriv_[num_nodes_ * num_int_pts_ * dim_];
};

}  // namespace nimble

#endif  // NIMBLE_ELEMENT_H
