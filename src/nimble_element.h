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
      double                                         density,
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

  static NIMBLE_INLINE_FUNCTION double
  Invert3x3(const double mat[][3], double inv[][3])
  {
    double minor0 = mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1];
    double minor1 = mat[1][0] * mat[2][2] - mat[1][2] * mat[2][0];
    double minor2 = mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0];
    double minor3 = mat[0][1] * mat[2][2] - mat[0][2] * mat[2][1];
    double minor4 = mat[0][0] * mat[2][2] - mat[2][0] * mat[0][2];
    double minor5 = mat[0][0] * mat[2][1] - mat[0][1] * mat[2][0];
    double minor6 = mat[0][1] * mat[1][2] - mat[0][2] * mat[1][1];
    double minor7 = mat[0][0] * mat[1][2] - mat[0][2] * mat[1][0];
    double minor8 = mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0];
    double det    = mat[0][0] * minor0 - mat[0][1] * minor1 + mat[0][2] * minor2;

    if (det <= 0.0) {
#ifdef NIMBLE_HAVE_KOKKOS
      if (det == 0.0)
        printf("\n**** Error in HexElement::Invert3x3(), singular matrix.\n");
      else
        printf(
            "\n**** Error in HexElement::Invert3x3(), negative determinant "
            "(%e)\n",
            det);
#else
      NIMBLE_ASSERT(det > 0.0, "\n**** Error in HexElement::Invert3x3(), singular matrix.\n");
#endif
    }

    inv[0][0] = minor0 / det;
    inv[0][1] = -1.0 * minor3 / det;
    inv[0][2] = minor6 / det;
    inv[1][0] = -1.0 * minor1 / det;
    inv[1][1] = minor4 / det;
    inv[1][2] = -1.0 * minor7 / det;
    inv[2][0] = minor2 / det;
    inv[2][1] = -1.0 * minor5 / det;
    inv[2][2] = minor8 / det;

    return det;
  }

  static NIMBLE_INLINE_FUNCTION void
  LU_Decompose(double mat[][3], int index[]);

  static NIMBLE_INLINE_FUNCTION void
  LU_Solve(const double a[][3], const int index[3], double b[3]);

  static NIMBLE_INLINE_FUNCTION void
  LU_Invert(const double mat[][3], double inv[][3]);

  static NIMBLE_INLINE_FUNCTION double
  MatrixInverseCheckCorrectness(const double mat[][3], const double inv[][3]);

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
  template <class ViewT>
  void
  ComputeConsistentMass_impl(
      double      density,
      const ViewT node_reference_coords,
      double      consistent_mass_matrix[][num_nodes_]) const
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
      double                                         density,
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

 protected:

  template< class ConstViewT, class TensorViewT >
  NIMBLE_FUNCTION
  void
  ComputeDeformationGradients_impl(
      ConstViewT node_reference_coords, ConstViewT node_displacements,
      TensorViewT  deformation_gradients) const
  {
    double rc1, rc2, rc3, cc1, cc2, cc3, sfd1, sfd2, sfd3;
    double b_inv[][3]    = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
    double def_grad[][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};

    // Loop over the integration points
    for (int int_pt = 0; int_pt < 8; int_pt++) {
      // \sum_{i}^{N_{node}} x_{i} \frac{\partial N_{i} (\xi)}{\partial \xi}
      double a[][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};

      // \sum_{i}^{N_{node}} X_{i} \frac{\partial N_{i} (\xi)}{\partial \xi}
      double b[][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};

      // Sum over the number of nodes
      for (int j = 0; j < 8; j++) {
        rc1 = node_reference_coords(j, 0);
        rc2 = node_reference_coords(j, 1);
        rc3 = node_reference_coords(j, 2);

        cc1 = node_reference_coords(j, 0) + node_displacements(j, 0);
        cc2 = node_reference_coords(j, 1) + node_displacements(j, 1);
        cc3 = node_reference_coords(j, 2) + node_displacements(j, 2);

        sfd1 = shape_fcn_deriv_[24 * int_pt + 3 * j];
        sfd2 = shape_fcn_deriv_[24 * int_pt + 3 * j + 1];
        sfd3 = shape_fcn_deriv_[24 * int_pt + 3 * j + 2];

        a[0][0] += cc1 * sfd1;
        a[0][1] += cc1 * sfd2;
        a[0][2] += cc1 * sfd3;
        a[1][0] += cc2 * sfd1;
        a[1][1] += cc2 * sfd2;
        a[1][2] += cc2 * sfd3;
        a[2][0] += cc3 * sfd1;
        a[2][1] += cc3 * sfd2;
        a[2][2] += cc3 * sfd3;

        b[0][0] += rc1 * sfd1;
        b[0][1] += rc1 * sfd2;
        b[0][2] += rc1 * sfd3;
        b[1][0] += rc2 * sfd1;
        b[1][1] += rc2 * sfd2;
        b[1][2] += rc2 * sfd3;
        b[2][0] += rc3 * sfd1;
        b[2][1] += rc3 * sfd2;
        b[2][2] += rc3 * sfd3;
      }

      Invert3x3(b, b_inv);

      for (int j = 0; j < 3; j++) {
        for (int k = 0; k < 3; k++) {
          def_grad[j][k] = a[j][0] * b_inv[0][k] + a[j][1] * b_inv[1][k] + a[j][2] * b_inv[2][k];
        }
      }

      deformation_gradients(int_pt, K_F_XX_) = def_grad[0][0];
      deformation_gradients(int_pt, K_F_XY_) = def_grad[0][1];
      deformation_gradients(int_pt, K_F_XZ_) = def_grad[0][2];
      deformation_gradients(int_pt, K_F_YX_) = def_grad[1][0];
      deformation_gradients(int_pt, K_F_YY_) = def_grad[1][1];
      deformation_gradients(int_pt, K_F_YZ_) = def_grad[1][2];
      deformation_gradients(int_pt, K_F_ZX_) = def_grad[2][0];
      deformation_gradients(int_pt, K_F_ZY_) = def_grad[2][1];
      deformation_gradients(int_pt, K_F_ZZ_) = def_grad[2][2];
    }
  }

 public:

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

 protected:

  template <class ConstViewT, class ViewTensorT, class ViewT>
  void
  ComputeNodalForces_impl(
      ConstViewT  node_reference_coords,
      ConstViewT  node_displacements,
      ViewTensorT int_pt_stresses,
      ViewT       node_forces) const
  {
    double cc1, cc2, cc3, sfd1, sfd2, sfd3, dN_dx1, dN_dx2, dN_dx3, f1, f2, f3;
    double jac_det;
    double a_inv[][dim_] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
    double force[][dim_] = {
        {0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0}};
    constexpr int dim_nodes = dim_ * num_nodes_;

    for (int int_pt = 0; int_pt < num_int_pts_; int_pt++) {
      // \sum_{i}^{N_{node}} x_{i} \frac{\partial N_{i} (\xi)}{\partial \xi}
      double a[][dim_] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};

      for (int n = 0; n < num_nodes_; n++) {
        cc1  = node_reference_coords(n, 0) + node_displacements(n, 0);
        cc2  = node_reference_coords(n, 1) + node_displacements(n, 1);
        cc3  = node_reference_coords(n, 2) + node_displacements(n, 2);
        sfd1 = shape_fcn_deriv_[dim_nodes * int_pt + dim_ * n];
        sfd2 = shape_fcn_deriv_[dim_nodes * int_pt + dim_ * n + 1];
        sfd3 = shape_fcn_deriv_[dim_nodes * int_pt + dim_ * n + 2];
        a[0][0] += cc1 * sfd1;
        a[0][1] += cc1 * sfd2;
        a[0][2] += cc1 * sfd3;
        a[1][0] += cc2 * sfd1;
        a[1][1] += cc2 * sfd2;
        a[1][2] += cc2 * sfd3;
        a[2][0] += cc3 * sfd1;
        a[2][1] += cc3 * sfd2;
        a[2][2] += cc3 * sfd3;
      }

      jac_det = Invert3x3(a, a_inv);

      for (int node = 0; node < num_nodes_; node++) {
        dN_dx1 = shape_fcn_deriv_[int_pt * dim_nodes + dim_ * node] * a_inv[0][0] +
                 shape_fcn_deriv_[int_pt * dim_nodes + dim_ * node + 1] * a_inv[1][0] +
                 shape_fcn_deriv_[int_pt * dim_nodes + dim_ * node + 2] * a_inv[2][0];

        dN_dx2 = shape_fcn_deriv_[int_pt * dim_nodes + dim_ * node] * a_inv[0][1] +
                 shape_fcn_deriv_[int_pt * dim_nodes + dim_ * node + 1] * a_inv[1][1] +
                 shape_fcn_deriv_[int_pt * dim_nodes + dim_ * node + 2] * a_inv[2][1];

        dN_dx3 = shape_fcn_deriv_[int_pt * dim_nodes + dim_ * node] * a_inv[0][2] +
                 shape_fcn_deriv_[int_pt * dim_nodes + dim_ * node + 1] * a_inv[1][2] +
                 shape_fcn_deriv_[int_pt * dim_nodes + dim_ * node + 2] * a_inv[2][2];

        f1 = dN_dx1 * int_pt_stresses(int_pt, K_S_XX_) + dN_dx2 * int_pt_stresses(int_pt, K_S_YX_) +
             dN_dx3 * int_pt_stresses(int_pt, K_S_ZX_);

        f2 = dN_dx1 * int_pt_stresses(int_pt, K_S_XY_) + dN_dx2 * int_pt_stresses(int_pt, K_S_YY_) +
             dN_dx3 * int_pt_stresses(int_pt, K_S_ZY_);

        f3 = dN_dx1 * int_pt_stresses(int_pt, K_S_XZ_) + dN_dx2 * int_pt_stresses(int_pt, K_S_YZ_) +
             dN_dx3 * int_pt_stresses(int_pt, K_S_ZZ_);

        f1 *= jac_det * int_wts_[int_pt];
        f2 *= jac_det * int_wts_[int_pt];
        f3 *= jac_det * int_wts_[int_pt];

        // profiling showed that calling the -= operator directly on the kokkos
        // view is expensive
        force[node][0] -= f1;
        force[node][1] -= f2;
        force[node][2] -= f3;
      }
    }

    for (int node = 0; node < num_nodes_; node++) {
      node_forces(node, 0) = force[node][0];
      node_forces(node, 1) = force[node][1];
      node_forces(node, 2) = force[node][2];
    }
  }

 public:
  void
  ComputeNodalForces(const double* node_current_coords, const double* int_pt_stresses, double* node_forces) override;

#ifdef NIMBLE_HAVE_KOKKOS
  NIMBLE_FUNCTION
  void
  ComputeNodalForces(
      nimble_kokkos::DeviceVectorNodeGatheredSubView node_reference_coords,
      nimble_kokkos::DeviceVectorNodeGatheredSubView node_displacements,
      nimble_kokkos::DeviceSymTensorIntPtSubView     int_pt_stresses,
      nimble_kokkos::DeviceVectorNodeGatheredSubView node_forces) const override
  {
    ComputeNodalForces_impl(node_reference_coords, node_displacements, int_pt_stresses, node_forces);
  }
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
