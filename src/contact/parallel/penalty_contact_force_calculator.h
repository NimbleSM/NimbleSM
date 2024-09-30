//
// Created by pierrelp on 9/26/24.
//

#ifndef NIMBLE_PENALTY_CONTACT_FORCE_CALCULATOR_H
#define NIMBLE_PENALTY_CONTACT_FORCE_CALCULATOR_H

#include "nimble_contact_force_calculator.h"

namespace nimble {

  class PenaltyContactForceCalculator : public ContactForceCalculator {
  public:
    void ComputeContactForce(double penalty, double gap, const double normal[3], double mass1, std::array<double, 3> a1, std::array<double, 3> a2, std::array<double, 3> N2, std::array<double, 3>& contact_force) const override {

      std::array<double, 3> predict_force = {0, 0, 0};
      std::array<double, 3> correct_force = {0, 0, 0};

      for (int i = 0; i < 3; ++i) {
        predict_force[i] = -penalty * gap * normal[i];
      }

      std::array<double, 3> cross_product;
      cross_product[0] = N2[1] * a2[2] - N2[2] * a2[1];
      cross_product[1] = N2[2] * a2[0] - N2[0] * a2[2];
      cross_product[2] = N2[0] * a2[1] - N2[1] * a2[0];

      for (int i = 0; i < 3; ++i) {
        double correction_factor = std::max(0.0, (cross_product[i] - a1[i]) * normal[i]);
        correct_force[i] = -mass1 * correction_factor * normal[i];
      }

      for (int i = 0; i < 3; ++i) {
        contact_force[i] = predict_force[i] + correct_force[i];
      }
    }
  };

} // namespace nimble

#endif //NIMBLE_PENALTY_CONTACT_FORCE_CALCULATOR_H
