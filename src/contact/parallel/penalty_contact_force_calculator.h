//
// Created by pierrelp on 9/26/24.
//

#ifndef NIMBLE_PENALTY_CONTACT_FORCE_CALCULATOR_H
#define NIMBLE_PENALTY_CONTACT_FORCE_CALCULATOR_H

#include "nimble_utils.h"
#include "nimble_contact_force_calculator.h"

namespace nimble {

  class PenaltyContactForceCalculator : public ContactForceCalculator {
  public:
    void ComputeContactForce(double penalty, double gap, const double normal[3], double mass1, double a1[3], double a2[3], double N2[3], std::array<double, 3>& contact_force) const override {

      double predict_force[3] = {0, 0, 0};
      double correct_force[3] = {0, 0, 0};

      for (int i = 0; i < 3; ++i) {
        predict_force[i] = -penalty * gap * normal[i];
      }

      double cross_product[3];
      CrossProduct(N2, a2, cross_product);

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
