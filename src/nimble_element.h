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
  Element() = default;

  NIMBLE_FUNCTION
  virtual ~Element() = default;

  virtual int
  Dim() const = 0;

  virtual int
  NumNodesPerElement() const = 0;

  virtual int
  NumIntegrationPointsPerElement() const = 0;

  virtual void
  ComputeLumpedMass(double density, const double* node_reference_coords, double* lumped_mass) const = 0;

#ifdef NIMBLE_HAVE_KOKKOS
  NIMBLE_FUNCTION
  virtual void
  ComputeLumpedMass(
      double                                   density,
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
  Invert3x3(const double mat[][3], double inv[][3]) const;

  static NIMBLE_INLINE_FUNCTION
  void
  LU_Decompose(double mat[][3], int index[]) ;

  static NIMBLE_INLINE_FUNCTION
  void
  LU_Solve(const double a[][3], const int index[3], double b[3]) ;

  static NIMBLE_INLINE_FUNCTION
  void
  LU_Invert(const double mat[][3], double inv[][3]) ;

  static NIMBLE_INLINE_FUNCTION
  double
  MatrixInverseCheckCorrectness(const double mat[][3], const double inv[][3]) ;

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
 private:
  static constexpr int dim_         = 3;
  static constexpr int num_nodes_   = 8;
  static constexpr int num_int_pts_ = 8;

 public:
  NIMBLE_FUNCTION
  HexElement();

  NIMBLE_FUNCTION
  ~HexElement() override = default;

  int
  Dim() const override
  {
    return dim_;
  }

  int
  NumNodesPerElement() const override
  {
    return num_nodes_;
  }

  int
  NumIntegrationPointsPerElement() const override
  {
    return num_int_pts_;
  }

 protected:

  template< class ViewT >
  void ComputeConsistentMass_impl(double density, const ViewT node_reference_coords,
                                  double consistent_mass_matrix[][num_nodes_]) const
  {
    double jac_det[num_int_pts_];
    double rc1, rc2, rc3, sfd1, sfd2, sfd3;
    for (int int_pt = 0; int_pt < num_int_pts_; int_pt++) {
      // \sum_{i}^{N_{node}} x_{i} \frac{\partial N_{i} (\xi)}{\partial \xi}
      double a[][dim_]     = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
      double a_inv[][dim_] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};

      for (int n = 0; n < num_nodes_; n++) {
        rc1  = node_reference_coords(n, 0);
        rc2  = node_reference_coords(n, 1);
        rc3  = node_reference_coords(n, 2);
        sfd1 = shape_fcn_deriv_[24 * int_pt + dim_ * n];
        sfd2 = shape_fcn_deriv_[24 * int_pt + dim_ * n + 1];
        sfd3 = shape_fcn_deriv_[24 * int_pt + dim_ * n + 2];
        a[0][0] += rc1 * sfd1;
        a[0][1] += rc1 * sfd2;
        a[0][2] += rc1 * sfd3;
        a[1][0] += rc2 * sfd1;
        a[1][1] += rc2 * sfd2;
        a[1][2] += rc2 * sfd3;
        a[2][0] += rc3 * sfd1;
        a[2][1] += rc3 * sfd2;
        a[2][2] += rc3 * sfd3;
      }
      jac_det[int_pt] = Invert3x3(a, a_inv);
    }

    for (int i = 0; i < num_nodes_; i++) {
      for (int j = 0; j < num_nodes_; j++) {
        consistent_mass_matrix[i][j] = 0.0;
        for (int int_pt = 0; int_pt < num_int_pts_; int_pt++) {
          consistent_mass_matrix[i][j] += int_wts_[int_pt] * density * shape_fcn_vals_[int_pt * num_nodes_ + i] *
              shape_fcn_vals_[int_pt * num_nodes_ + j] * jac_det[int_pt];
        }
      }
    }
  }

 public:

  void
  ComputeLumpedMass(double density, const double* node_reference_coords, double* lumped_mass) const override;

#ifdef NIMBLE_HAVE_KOKKOS
  NIMBLE_FUNCTION
  void
  ComputeLumpedMass(
      double                                   density,
      nimble_kokkos::DeviceVectorNodeGatheredSubView node_reference_coords,
      nimble_kokkos::DeviceScalarNodeGatheredSubView lumped_mass) const override;
#endif

  double
  ComputeCharacteristicLength(const double* node_coords) override;

  void
  ComputeVolumeAverage(
      const double* node_current_coords,
      int           num_quantities,
      const double* int_pt_quantities,
      double&       volume,
      double*       volume_averaged_quantity) const override;

#ifdef NIMBLE_HAVE_KOKKOS

  NIMBLE_FUNCTION
  void
  ComputeVolume(
      nimble_kokkos::DeviceVectorNodeGatheredSubView node_reference_coords,
      nimble_kokkos::DeviceVectorNodeGatheredSubView node_displacements,
      nimble_kokkos::DeviceScalarElemSingleEntryView elem_volume) const override;

  NIMBLE_FUNCTION
  void
  ComputeVolumeAverageSymTensor(
      nimble_kokkos::DeviceVectorNodeGatheredSubView    node_reference_coords,
      nimble_kokkos::DeviceVectorNodeGatheredSubView    node_displacements,
      nimble_kokkos::DeviceSymTensorIntPtSubView        int_pt_quantities,
      nimble_kokkos::DeviceSymTensorElemSingleEntryView vol_ave_quantity) const override;

  NIMBLE_FUNCTION
  void
  ComputeVolumeAverageFullTensor(
      nimble_kokkos::DeviceVectorNodeGatheredSubView     node_reference_coords,
      nimble_kokkos::DeviceVectorNodeGatheredSubView     node_displacements,
      nimble_kokkos::DeviceFullTensorIntPtSubView        int_pt_quantities,
      nimble_kokkos::DeviceFullTensorElemSingleEntryView vol_ave_quantity) const override;
#endif

  void
  ComputeDeformationGradients(
      const double* node_reference_coords,
      const double* node_current_coords,
      double*       deformation_gradients) const override;

#ifdef NIMBLE_HAVE_KOKKOS
  NIMBLE_FUNCTION
  void
  ComputeDeformationGradients(
      nimble_kokkos::DeviceVectorNodeGatheredSubView node_reference_coords,
      nimble_kokkos::DeviceVectorNodeGatheredSubView node_displacements,
      nimble_kokkos::DeviceFullTensorIntPtSubView    deformation_gradients) const override;
#endif

  void
  ComputeTangent(const double* node_current_coords, const double* material_tangent, double* element_tangent) override;

  void
  ComputeNodalForces(const double* node_current_coords, const double* int_pt_stresses, double* node_forces) override;

#ifdef NIMBLE_HAVE_KOKKOS
  NIMBLE_FUNCTION
  void
  ComputeNodalForces(
      nimble_kokkos::DeviceVectorNodeGatheredSubView node_reference_coords,
      nimble_kokkos::DeviceVectorNodeGatheredSubView node_displacements,
      nimble_kokkos::DeviceSymTensorIntPtSubView     int_pt_stresses,
      nimble_kokkos::DeviceVectorNodeGatheredSubView node_forces) const override;
#endif

 protected:
  NIMBLE_FUNCTION
  void
  ShapeFunctionValues(const double* natural_coords, double* shape_function_values);

  NIMBLE_FUNCTION
  void
  ShapeFunctionDerivatives(const double* natural_coords, double* shape_function_derivatives);

 private:

  double int_pts_[num_int_pts_ * dim_]{};
  double int_wts_[num_int_pts_]{};
  double shape_fcn_vals_[num_nodes_ * num_int_pts_]{};
  double shape_fcn_deriv_[num_nodes_ * num_int_pts_ * dim_]{};
};

}  // namespace nimble

#endif  // NIMBLE_ELEMENT_H
