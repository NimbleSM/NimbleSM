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

#include "nimble_element.h"
#include "nimble_utils.h"
#include <limits>

namespace nimble {

  double Element::Invert3x3(double mat[][3],
                            double inv[][3]) const {

    double minor0 =  mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1];
    double minor1 =  mat[1][0] * mat[2][2] - mat[1][2] * mat[2][0];
    double minor2 =  mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0];
    double minor3 =  mat[0][1] * mat[2][2] - mat[0][2] * mat[2][1];
    double minor4 =  mat[0][0] * mat[2][2] - mat[2][0] * mat[0][2];
    double minor5 =  mat[0][0] * mat[2][1] - mat[0][1] * mat[2][0];
    double minor6 =  mat[0][1] * mat[1][2] - mat[0][2] * mat[1][1];
    double minor7 =  mat[0][0] * mat[1][2] - mat[0][2] * mat[1][0];
    double minor8 =  mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0];
    double det = mat[0][0] * minor0 - mat[0][1] * minor1 + mat[0][2] * minor2;

    if (det <= 0.0) {
      printf("\n**** Error in HexElement::Invert3x3(), singular matrix.\n");
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

  HexElement::HexElement() {

    // 1/sqrt(3)
    const double val = 0.577350269189626 ;
    // integration point 1
    int_pts_[0]  = -val;
    int_pts_[1]  = -val;
    int_pts_[2]  = -val;
    // integration point 2
    int_pts_[3]  =  val;
    int_pts_[4]  = -val;
    int_pts_[5]  = -val;
    // integration point 3
    int_pts_[6]  =  val;
    int_pts_[7]  =  val;
    int_pts_[8]  = -val;
    // integration point 4
    int_pts_[9]  = -val;
    int_pts_[10] =  val;
    int_pts_[11] = -val;
    // integration point 5
    int_pts_[12] = -val;
    int_pts_[13] = -val;
    int_pts_[14] =  val;
    // integration point 6
    int_pts_[15] =  val;
    int_pts_[16] = -val;
    int_pts_[17] =  val;
    // integration point 7
    int_pts_[18] =  val;
    int_pts_[19] =  val;
    int_pts_[20] =  val;
    // integration point 8
    int_pts_[21] = -val;
    int_pts_[22] =  val;
    int_pts_[23] =  val;

    for (int i=0 ; i<num_int_pts_ ; i++) {
      int_wts_[i] = 1.0;
    }

    ShapeFunctionValues(int_pts_, shape_fcn_vals_);

    ShapeFunctionDerivatives(int_pts_, shape_fcn_deriv_);
  }

  void HexElement::ShapeFunctionValues(const double* natural_coords,
                                       double* shape_function_values) {
    double r, s, t;
    double c = 1.0/8.0;

    // Loop over the integration points
    for (int i=0 ; i<8 ; i++) {

      // Natural coordinates of this integration point
      r = natural_coords[3*i];
      s = natural_coords[3*i + 1];
      t = natural_coords[3*i + 2];

      // Value of each of the eight shape functions at this integration point
      shape_function_values[8*i]     = c * (1.0 - r) * (1.0 - s) * (1.0 - t);
      shape_function_values[8*i + 1] = c * (1.0 + r) * (1.0 - s) * (1.0 - t);
      shape_function_values[8*i + 2] = c * (1.0 + r) * (1.0 + s) * (1.0 - t);
      shape_function_values[8*i + 3] = c * (1.0 - r) * (1.0 + s) * (1.0 - t);
      shape_function_values[8*i + 4] = c * (1.0 - r) * (1.0 - s) * (1.0 + t);
      shape_function_values[8*i + 5] = c * (1.0 + r) * (1.0 - s) * (1.0 + t);
      shape_function_values[8*i + 6] = c * (1.0 + r) * (1.0 + s) * (1.0 + t);
      shape_function_values[8*i + 7] = c * (1.0 - r) * (1.0 + s) * (1.0 + t);
    }
  }

  void HexElement::ShapeFunctionDerivatives(const double* natural_coords,
                                            double* shape_function_derivatives) {
    double r, s, t;
    double c = 1.0/8.0;

    // Loop over the integration points
    for (int int_pt=0 ; int_pt<8 ; int_pt++) {

      // Natural coordinates of this integration point
      r = natural_coords[3*int_pt];
      s = natural_coords[3*int_pt + 1];
      t = natural_coords[3*int_pt + 2];

      // Derivative of each of the eight shape functions w.r.t. the natural coordinates at this integration point
      // shape function 1
      shape_function_derivatives[24*int_pt]      = -c * (1.0 - s) * (1.0 - t);
      shape_function_derivatives[24*int_pt + 1]  = -c * (1.0 - r) * (1.0 - t);
      shape_function_derivatives[24*int_pt + 2]  = -c * (1.0 - r) * (1.0 - s);
      // shape function 2
      shape_function_derivatives[24*int_pt + 3]  =  c * (1.0 - s) * (1.0 - t);
      shape_function_derivatives[24*int_pt + 4]  = -c * (1.0 + r) * (1.0 - t);
      shape_function_derivatives[24*int_pt + 5]  = -c * (1.0 + r) * (1.0 - s);
      // shape function 3
      shape_function_derivatives[24*int_pt + 6]  =  c * (1.0 + s) * (1.0 - t);
      shape_function_derivatives[24*int_pt + 7]  =  c * (1.0 + r) * (1.0 - t);
      shape_function_derivatives[24*int_pt + 8]  = -c * (1.0 + r) * (1.0 + s);
      // shape function 4
      shape_function_derivatives[24*int_pt + 9]  = -c * (1.0 + s) * (1.0 - t);
      shape_function_derivatives[24*int_pt + 10] =  c * (1.0 - r) * (1.0 - t);
      shape_function_derivatives[24*int_pt + 11] = -c * (1.0 - r) * (1.0 + s);
      // shape function 5
      shape_function_derivatives[24*int_pt + 12] = -c * (1.0 - s) * (1.0 + t);
      shape_function_derivatives[24*int_pt + 13] = -c * (1.0 - r) * (1.0 + t);
      shape_function_derivatives[24*int_pt + 14] =  c * (1.0 - r) * (1.0 - s);
      // shape function 6
      shape_function_derivatives[24*int_pt + 15] =  c * (1.0 - s) * (1.0 + t);
      shape_function_derivatives[24*int_pt + 16] = -c * (1.0 + r) * (1.0 + t);
      shape_function_derivatives[24*int_pt + 17] =  c * (1.0 + r) * (1.0 - s);
      // shape function 7
      shape_function_derivatives[24*int_pt + 18] =  c * (1.0 + s) * (1.0 + t);
      shape_function_derivatives[24*int_pt + 19] =  c * (1.0 + r) * (1.0 + t);
      shape_function_derivatives[24*int_pt + 20] =  c * (1.0 + r) * (1.0 + s);
      // shape function 8
      shape_function_derivatives[24*int_pt + 21] = -c * (1.0 + s) * (1.0 + t);
      shape_function_derivatives[24*int_pt + 22] =  c * (1.0 - r) * (1.0 + t);
      shape_function_derivatives[24*int_pt + 23] =  c * (1.0 - r) * (1.0 + s);
    }
  }

  void HexElement::ComputeLumpedMass(const double  density,
                                     const double* node_reference_coords,
                                     double* lumped_mass) const {

    double consistent_mass_matrix[][24] = {{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                                           {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                                           {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                                           {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                                           {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                                           {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                                           {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                                           {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};

    double jac_det[num_int_pts_];
    double rc1, rc2, rc3, sfd1, sfd2, sfd3;
    for (int int_pt=0 ; int_pt<num_int_pts_ ; int_pt++) {

      // \sum_{i}^{N_{node}} x_{i} \frac{\partial N_{i} (\xi)}{\partial \xi}
      double a[][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
      double a_inv[][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};

      for (int n=0 ; n<num_nodes_ ; n++) {
        rc1 = node_reference_coords[3*n];
        rc2 = node_reference_coords[3*n + 1];
        rc3 = node_reference_coords[3*n + 2];
        sfd1 = shape_fcn_deriv_[24*int_pt + 3*n];
        sfd2 = shape_fcn_deriv_[24*int_pt + 3*n + 1];
        sfd3 = shape_fcn_deriv_[24*int_pt + 3*n + 2];
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

    for (int i=0 ; i<8 ; i++) {
      for (int j=0 ; j<8 ; j++) {
        consistent_mass_matrix[i][j] = 0.0;
        for (int int_pt=0 ; int_pt<num_int_pts_ ; int_pt++) {
          consistent_mass_matrix[i][j] += int_wts_[int_pt] * density * shape_fcn_vals_[int_pt*num_nodes_ + i] * shape_fcn_vals_[int_pt*num_nodes_ + j] * jac_det[int_pt];
        }
      }
    }

    for (int i=0 ; i<8 ; i++) {
      lumped_mass[i] = 0.0;
      for (int j=0 ; j<8 ; j++) {
        lumped_mass[i] += consistent_mass_matrix[i][j];
      }
    }
  }

#ifdef NIMBLE_HAVE_KOKKOS
  void HexElement::ComputeLumpedMass(const double density,
                                     nimble_kokkos::DeviceVectorGatheredSubView node_reference_coords,
                                     nimble_kokkos::DeviceScalarGatheredSubView lumped_mass) const {

    double consistent_mass_matrix[][24] = {{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                                           {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                                           {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                                           {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                                           {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                                           {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                                           {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                                           {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};

    double jac_det[num_int_pts_];
    double rc1, rc2, rc3, sfd1, sfd2, sfd3;
    for (int int_pt=0 ; int_pt<num_int_pts_ ; int_pt++) {

      // \sum_{i}^{N_{node}} x_{i} \frac{\partial N_{i} (\xi)}{\partial \xi}
      double a[][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
      double a_inv[][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};

      for (int n=0 ; n<num_nodes_ ; n++) {
        rc1 = node_reference_coords(n,0);
        rc2 = node_reference_coords(n,1);
        rc3 = node_reference_coords(n,2);
        sfd1 = shape_fcn_deriv_[24*int_pt + 3*n];
        sfd2 = shape_fcn_deriv_[24*int_pt + 3*n + 1];
        sfd3 = shape_fcn_deriv_[24*int_pt + 3*n + 2];
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

    for (int i=0 ; i<8 ; i++) {
      for (int j=0 ; j<8 ; j++) {
        consistent_mass_matrix[i][j] = 0.0;
        for (int int_pt=0 ; int_pt<num_int_pts_ ; int_pt++) {
          consistent_mass_matrix[i][j] += int_wts_[int_pt] * density * shape_fcn_vals_[int_pt*num_nodes_ + i] * shape_fcn_vals_[int_pt*num_nodes_ + j] * jac_det[int_pt];
        }
      }
    }

    for (int i=0 ; i<8 ; i++) {
      lumped_mass(i) = 0.0;
      for (int j=0 ; j<8 ; j++) {
        lumped_mass(i) += consistent_mass_matrix[i][j];
      }
    }
  }
#endif

  double HexElement::ComputeCharacteristicLength(const double* node_coords) {

    // TODO Implement a better algorithm for finding the minimum
    //      length across the element.

    double characteristic_length, distance_squared, min_distance_squared, nx, ny, nz, mx, my, mz;
    double x_min, x_max, y_min, y_max, z_min, z_max;

    x_max = y_max = z_max = 0.0;
    min_distance_squared = x_min = y_min = z_min = std::numeric_limits<double>::max();

    for (int n=0 ; n<num_nodes_ ; n++) {
      nx = node_coords[3*n];
      ny = node_coords[3*n + 1];
      nz = node_coords[3*n + 2];
      if (nx < x_min) {
        x_min = nx;
      }
      if (nx > x_max) {
        x_max = nx;
      }
      if (ny < y_min) {
        y_min = ny;
      }
      if (ny > y_max) {
        y_max = ny;
      }
      if (nz < z_min) {
        z_min = nz;
      }
      if (nz > z_max) {
        z_max = nz;
      }
      for (int m=n+1 ; m<num_nodes_ ; m++) {
        mx = node_coords[3*m];
        my = node_coords[3*m + 1];
        mz = node_coords[3*m + 2];
        distance_squared = (nx - mx)*(nx - mx) + (ny - my)*(ny - my) + (nz - mz)*(nz - mz);
        if (distance_squared < min_distance_squared) {
          min_distance_squared = distance_squared;
        }
      }
    }
    characteristic_length = std::sqrt(min_distance_squared);

    double min_box_length = x_max - x_min;
    if (y_max - y_min < min_box_length) {
      min_box_length = y_max - y_min;
    }
    if (z_max - z_min < min_box_length) {
      min_box_length = z_max - z_min;
    }

    if (min_box_length < characteristic_length) {
      characteristic_length = min_box_length;
    }

    return characteristic_length;
  }

  void HexElement::ComputeVolumeAverage(const double* node_reference_coords,
                                        const double* node_current_coords,
                                        int num_quantities,
                                        const double* int_pt_quantities,
                                        double& volume,
                                        double* volume_averaged_quantities) {

    double cc1, cc2, cc3, sfd1, sfd2, sfd3;
    double jac_det;
    double a_inv[][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};

    volume = 0.0;
    for (int i_quantity = 0 ; i_quantity < num_quantities ; i_quantity++) {
      volume_averaged_quantities[i_quantity] = 0.0;
    }

    for (int int_pt=0 ; int_pt<num_int_pts_ ; int_pt++) {

      // \sum_{i}^{N_{node}} x_{i} \frac{\partial N_{i} (\xi)}{\partial \xi}
      double a[][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};

      for (int n=0 ; n<num_nodes_ ; n++) {
        cc1 = node_current_coords[3*n];
        cc2 = node_current_coords[3*n + 1];
        cc3 = node_current_coords[3*n + 2];
        sfd1 = shape_fcn_deriv_[24*int_pt + 3*n];
        sfd2 = shape_fcn_deriv_[24*int_pt + 3*n + 1];
        sfd3 = shape_fcn_deriv_[24*int_pt + 3*n + 2];
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
      volume += jac_det;
      for (int i_quantity = 0 ; i_quantity < num_quantities ; i_quantity++) {
        volume_averaged_quantities[i_quantity] += int_pt_quantities[int_pt*num_quantities + i_quantity]*int_wts_[int_pt]*jac_det;
      }
    }

    for (int i_quantity = 0 ; i_quantity < num_quantities ; i_quantity++) {
      volume_averaged_quantities[i_quantity] /= volume;
    }
  }

  void HexElement::ComputeDeformationGradients(const double* node_reference_coords,
                                               const double* node_current_coords,
                                               double* deformation_gradients) const {

    double rc1, rc2, rc3, cc1, cc2, cc3, sfd1, sfd2, sfd3;
    double b_inv[][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
    double def_grad[][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};

    // Loop over the integration points
    for (int int_pt=0 ; int_pt<8 ; int_pt++) {

      // \sum_{i}^{N_{node}} x_{i} \frac{\partial N_{i} (\xi)}{\partial \xi}
      double a[][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};

      // \sum_{i}^{N_{node}} X_{i} \frac{\partial N_{i} (\xi)}{\partial \xi}
      double b[][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};

      // Sum over the number of nodes
      for (int j=0 ; j<8 ; j++) {

        rc1 = node_reference_coords[3*j];
        rc2 = node_reference_coords[3*j + 1];
        rc3 = node_reference_coords[3*j + 2];

        cc1 = node_current_coords[3*j];
        cc2 = node_current_coords[3*j + 1];
        cc3 = node_current_coords[3*j + 2];

        sfd1 = shape_fcn_deriv_[24*int_pt + 3*j];
        sfd2 = shape_fcn_deriv_[24*int_pt + 3*j + 1];
        sfd3 = shape_fcn_deriv_[24*int_pt + 3*j + 2];

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

      for (int j=0 ; j<3 ; j++) {
        for (int k=0 ; k<3 ; k++) {
          def_grad[j][k] = a[j][0] * b_inv[0][k] + a[j][1] * b_inv[1][k] + a[j][2] * b_inv[2][k];
        }
      }

      deformation_gradients[9*int_pt + K_F_XX] = def_grad[0][0];
      deformation_gradients[9*int_pt + K_F_XY] = def_grad[0][1];
      deformation_gradients[9*int_pt + K_F_XZ] = def_grad[0][2];
      deformation_gradients[9*int_pt + K_F_YX] = def_grad[1][0];
      deformation_gradients[9*int_pt + K_F_YY] = def_grad[1][1];
      deformation_gradients[9*int_pt + K_F_YZ] = def_grad[1][2];
      deformation_gradients[9*int_pt + K_F_ZX] = def_grad[2][0];
      deformation_gradients[9*int_pt + K_F_ZY] = def_grad[2][1];
      deformation_gradients[9*int_pt + K_F_ZZ] = def_grad[2][2];
    }
  }

#ifdef NIMBLE_HAVE_KOKKOS
  void HexElement::ComputeDeformationGradients(nimble_kokkos::DeviceVectorGatheredSubView node_reference_coords,
                                               nimble_kokkos::DeviceVectorGatheredSubView node_displacements,
                                               nimble_kokkos::DeviceFullTensorSubView deformation_gradients) const {

    double rc1, rc2, rc3, cc1, cc2, cc3, sfd1, sfd2, sfd3;
    double b_inv[][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
    double def_grad[][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};

    // Loop over the integration points
    for (int int_pt=0 ; int_pt<8 ; int_pt++) {

      // \sum_{i}^{N_{node}} x_{i} \frac{\partial N_{i} (\xi)}{\partial \xi}
      double a[][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};

      // \sum_{i}^{N_{node}} X_{i} \frac{\partial N_{i} (\xi)}{\partial \xi}
      double b[][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};

      // Sum over the number of nodes
      for (int j=0 ; j<8 ; j++) {

        rc1 = node_reference_coords(j,0);
        rc2 = node_reference_coords(j,1);
        rc3 = node_reference_coords(j,2);

        cc1 = node_reference_coords(j,0) + node_displacements(j,0);
        cc2 = node_reference_coords(j,1) + node_displacements(j,1);
        cc3 = node_reference_coords(j,2) + node_displacements(j,2);

        sfd1 = shape_fcn_deriv_[24*int_pt + 3*j];
        sfd2 = shape_fcn_deriv_[24*int_pt + 3*j + 1];
        sfd3 = shape_fcn_deriv_[24*int_pt + 3*j + 2];

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

      for (int j=0 ; j<3 ; j++) {
        for (int k=0 ; k<3 ; k++) {
          def_grad[j][k] = a[j][0] * b_inv[0][k] + a[j][1] * b_inv[1][k] + a[j][2] * b_inv[2][k];
        }
      }

      deformation_gradients(int_pt,K_F_XX_) = def_grad[0][0];
      deformation_gradients(int_pt,K_F_XY_) = def_grad[0][1];
      deformation_gradients(int_pt,K_F_XZ_) = def_grad[0][2];
      deformation_gradients(int_pt,K_F_YX_) = def_grad[1][0];
      deformation_gradients(int_pt,K_F_YY_) = def_grad[1][1];
      deformation_gradients(int_pt,K_F_YZ_) = def_grad[1][2];
      deformation_gradients(int_pt,K_F_ZX_) = def_grad[2][0];
      deformation_gradients(int_pt,K_F_ZY_) = def_grad[2][1];
      deformation_gradients(int_pt,K_F_ZZ_) = def_grad[2][2];
    }
  }
#endif

  void HexElement::ComputeTangent(const double* node_current_coords,
                                  const double* material_tangent,
                                  double* element_tangent) {

    double jac_det;
    double cc1, cc2, cc3, sfd1, sfd2, sfd3;

    for (int i=0 ; i<24 ; i++) {
      for (int j=0 ; j<24 ; j++) {
        element_tangent[i*24 + j] = 0.0;
      }
    }

    // \mathbf{K}_{elem} = \int_{\Omega_{o}} \mathbf{B}^{T}_{0} \mathbf{C}^{SE} \mathbf{B}_{o} \Omega_{o}

    //      dN1/dX    0       0     dN2/dX    0       0     ...  dN8/dX    0       0
    //        0     dN1/dY    0       0     dN2/dY    0     ...    0     dN8/dY    0
    //        0       0     dN1/dZ    0       0     dN2/dZ  ...    0       0     dN8/dZ
    // B =  dN1/dY  dN1/dX    0     dN2/dY  dN2/dX    0     ...  dN8/dY  dN8/dX    0
    //        0     dN1/dZ  dN1/dY    0     dN2/dZ  dN2/dY  ...    0     dN8/dZ  dN8/dY
    //      dN1/dZ    0     dN1/dX  dN2/dZ    0     dN2/dX  ...  dN8/dZ    0     dN8/dX

    double B[6][24];
    for (int i=0 ; i<6 ; i++) {
      for (int j=0 ; j<24 ; j++) {
        B[i][j] = 0.0;
      }
    }

    for (int int_pt=0 ; int_pt<num_int_pts_ ; int_pt++) {

      // \mathbf{a} = \sum_{i}^{N_{node}} x_{i} \frac{\partial N_{i} (\xi)}{\partial \xi}
      double a[][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
      double a_inv[][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
      double temp[6][24];
      for (int i=0 ; i<6 ; i++) {
        for (int j=0 ; j<24 ; j++) {
          temp[i][j] = 0.0;
        }
      }

      for (int n=0 ; n<num_nodes_ ; n++) {
        cc1 = node_current_coords[3*n];
        cc2 = node_current_coords[3*n + 1];
        cc3 = node_current_coords[3*n + 2];
        sfd1 = shape_fcn_deriv_[24*int_pt + 3*n];
        sfd2 = shape_fcn_deriv_[24*int_pt + 3*n + 1];
        sfd3 = shape_fcn_deriv_[24*int_pt + 3*n + 2];
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

      // derivatives of shape function with respect to current coordinate
      double dN_dcc1, dN_dcc2, dN_dcc3;
      for (int n=0 ; n<num_nodes_ ; n++) {
        sfd1 = shape_fcn_deriv_[24*int_pt + 3*n];
        sfd2 = shape_fcn_deriv_[24*int_pt + 3*n + 1];
        sfd3 = shape_fcn_deriv_[24*int_pt + 3*n + 2];
        dN_dcc1 = sfd1 * a_inv[0][0] + sfd2 * a_inv[1][0] + sfd3 * a_inv[2][0];
        dN_dcc2 = sfd1 * a_inv[0][1] + sfd2 * a_inv[1][1] + sfd3 * a_inv[2][1];
        dN_dcc3 = sfd1 * a_inv[0][2] + sfd2 * a_inv[1][2] + sfd3 * a_inv[2][2];
        B[0][3*n]     = dN_dcc1;
        B[1][3*n + 1] = dN_dcc2;
        B[2][3*n + 2] = dN_dcc3;
        B[3][3*n]     = dN_dcc2;
        B[3][3*n + 1] = dN_dcc1;
        B[4][3*n + 1] = dN_dcc3;
        B[4][3*n + 2] = dN_dcc2;
        B[5][3*n]     = dN_dcc3;
        B[5][3*n + 2] = dN_dcc1;
      }

      for (int i=0 ; i<6 ; i++) {
        for (int j=0 ; j<24 ; j++) {
          for (int k=0 ; k<6 ; k++) {
            temp[i][j] += material_tangent[36*int_pt + 6*i + k] * B[k][j];
          }
        }
      }

      for (int i=0 ; i<24 ; i++) {
        for (int j=0 ; j<24 ; j++) {
          for (int k=0 ; k<6 ; k++) {
            element_tangent[i*24 + j] += B[k][i] * temp[k][j] * int_wts_[int_pt] * jac_det;
          }
        }
      }
    }
  }

  void HexElement::ComputeNodalForces(const double* node_current_coords,
                                      const double* int_pt_stresses,
                                      double* node_forces) {

    double cc1, cc2, cc3, sfd1, sfd2, sfd3, dN_dx1, dN_dx2, dN_dx3, f1, f2, f3;
    double jac_det;
    double a_inv[][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
    int sym_tensor_size = 6;

    for (int n=0 ; n<num_nodes_; n++) {
      for (int i=0 ; i<dim_ ; i++) {
        node_forces[n*dim_ + i] = 0.0;
      }
    }

    for (int int_pt=0 ; int_pt<num_int_pts_ ; int_pt++) {

      // \sum_{i}^{N_{node}} x_{i} \frac{\partial N_{i} (\xi)}{\partial \xi}
      double a[][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};

      for (int n=0 ; n<num_nodes_ ; n++) {
        cc1 = node_current_coords[3*n];
        cc2 = node_current_coords[3*n + 1];
        cc3 = node_current_coords[3*n + 2];
        sfd1 = shape_fcn_deriv_[24*int_pt + 3*n];
        sfd2 = shape_fcn_deriv_[24*int_pt + 3*n + 1];
        sfd3 = shape_fcn_deriv_[24*int_pt + 3*n + 2];
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

      for (int node=0 ; node<num_nodes_ ; node++) {

        dN_dx1 =
          shape_fcn_deriv_[int_pt*24 + 3*node]     * a_inv[0][0] +
          shape_fcn_deriv_[int_pt*24 + 3*node + 1] * a_inv[1][0] +
          shape_fcn_deriv_[int_pt*24 + 3*node + 2] * a_inv[2][0];

        dN_dx2 =
          shape_fcn_deriv_[int_pt*24 + 3*node]     * a_inv[0][1] +
          shape_fcn_deriv_[int_pt*24 + 3*node + 1] * a_inv[1][1] +
          shape_fcn_deriv_[int_pt*24 + 3*node + 2] * a_inv[2][1];

        dN_dx3 =
          shape_fcn_deriv_[int_pt*24 + 3*node]     * a_inv[0][2] +
          shape_fcn_deriv_[int_pt*24 + 3*node + 1] * a_inv[1][2] +
          shape_fcn_deriv_[int_pt*24 + 3*node + 2] * a_inv[2][2];

        f1 =
          dN_dx1 * int_pt_stresses[sym_tensor_size*int_pt + K_S_XX] +
          dN_dx2 * int_pt_stresses[sym_tensor_size*int_pt + K_S_YX] +
          dN_dx3 * int_pt_stresses[sym_tensor_size*int_pt + K_S_ZX];

        f2 =
          dN_dx1 * int_pt_stresses[sym_tensor_size*int_pt + K_S_XY] +
          dN_dx2 * int_pt_stresses[sym_tensor_size*int_pt + K_S_YY] +
          dN_dx3 * int_pt_stresses[sym_tensor_size*int_pt + K_S_ZY];

        f3 =
          dN_dx1 * int_pt_stresses[sym_tensor_size*int_pt + K_S_XZ] +
          dN_dx2 * int_pt_stresses[sym_tensor_size*int_pt + K_S_YZ] +
          dN_dx3 * int_pt_stresses[sym_tensor_size*int_pt + K_S_ZZ];

        f1 *= jac_det * int_wts_[int_pt];
        f2 *= jac_det * int_wts_[int_pt];
        f3 *= jac_det * int_wts_[int_pt];

        node_forces[3*node    ] -= f1;
        node_forces[3*node + 1] -= f2;
        node_forces[3*node + 2] -= f3;
      }
    }
  }

#ifdef NIMBLE_HAVE_KOKKOS
  void HexElement::ComputeNodalForces(nimble_kokkos::DeviceVectorGatheredSubView node_reference_coords,
                                      nimble_kokkos::DeviceVectorGatheredSubView node_displacements,
                                      nimble_kokkos::DeviceSymTensorSubView int_pt_stresses,
                                      nimble_kokkos::DeviceVectorGatheredSubView node_forces) const {

    double cc1, cc2, cc3, sfd1, sfd2, sfd3, dN_dx1, dN_dx2, dN_dx3, f1, f2, f3;
    double jac_det;
    double a_inv[][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
    double force[][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};

    for (int int_pt=0 ; int_pt<num_int_pts_ ; int_pt++) {

      // \sum_{i}^{N_{node}} x_{i} \frac{\partial N_{i} (\xi)}{\partial \xi}
      double a[][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};

      for (int n=0 ; n<num_nodes_ ; n++) {
        cc1 = node_reference_coords(n,0) + node_displacements(n,0);
        cc2 = node_reference_coords(n,1) + node_displacements(n,1);
        cc3 = node_reference_coords(n,2) + node_displacements(n,2);
        sfd1 = shape_fcn_deriv_[24*int_pt + 3*n];
        sfd2 = shape_fcn_deriv_[24*int_pt + 3*n + 1];
        sfd3 = shape_fcn_deriv_[24*int_pt + 3*n + 2];
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

      for (int node=0 ; node<num_nodes_ ; node++) {

        dN_dx1 =
          shape_fcn_deriv_[int_pt*24 + 3*node]     * a_inv[0][0] +
          shape_fcn_deriv_[int_pt*24 + 3*node + 1] * a_inv[1][0] +
          shape_fcn_deriv_[int_pt*24 + 3*node + 2] * a_inv[2][0];

        dN_dx2 =
          shape_fcn_deriv_[int_pt*24 + 3*node]     * a_inv[0][1] +
          shape_fcn_deriv_[int_pt*24 + 3*node + 1] * a_inv[1][1] +
          shape_fcn_deriv_[int_pt*24 + 3*node + 2] * a_inv[2][1];

        dN_dx3 =
          shape_fcn_deriv_[int_pt*24 + 3*node]     * a_inv[0][2] +
          shape_fcn_deriv_[int_pt*24 + 3*node + 1] * a_inv[1][2] +
          shape_fcn_deriv_[int_pt*24 + 3*node + 2] * a_inv[2][2];

        f1 =
          dN_dx1 * int_pt_stresses(int_pt, K_S_XX_) +
          dN_dx2 * int_pt_stresses(int_pt, K_S_YX_) +
          dN_dx3 * int_pt_stresses(int_pt, K_S_ZX_);

        f2 =
          dN_dx1 * int_pt_stresses(int_pt, K_S_XY_) +
          dN_dx2 * int_pt_stresses(int_pt, K_S_YY_) +
          dN_dx3 * int_pt_stresses(int_pt, K_S_ZY_);

        f3 =
          dN_dx1 * int_pt_stresses(int_pt, K_S_XZ_) +
          dN_dx2 * int_pt_stresses(int_pt, K_S_YZ_) +
          dN_dx3 * int_pt_stresses(int_pt, K_S_ZZ_);

        f1 *= jac_det * int_wts_[int_pt];
        f2 *= jac_det * int_wts_[int_pt];
        f3 *= jac_det * int_wts_[int_pt];

        // profiling showed that calling the -= operator directly on the kokkos view is expensive
        force[node][0] -= f1;
        force[node][1] -= f2;
        force[node][2] -= f3;
      }
    }

    for (int node=0 ; node<num_nodes_ ; node++) {
      node_forces(node, 0) = force[node][0];
      node_forces(node, 1) = force[node][1];
      node_forces(node, 2) = force[node][2];
    }
  }
#endif

} // namespace nimble
