//
// Created by pierrelp on 9/26/24.
//

#ifndef NIMBLE_PENALTY_CONTACT_FORCE_CALCULATOR_H
#define NIMBLE_PENALTY_CONTACT_FORCE_CALCULATOR_H

#include "nimble_contact_force_calculator.h"

namespace nimble {

  class PenaltyContactForceCalculator : public ContactForceCalculator {
  public:
    void ComputeContactForce(double penalty, double gap, const double normal[3], std::array<double, 3>& contact_force) const override {
      for (int i = 0; i < 3; ++i) {
        contact_force[i] = -penalty * gap * normal[i];
      }
    }
  };

} // namespace nimble

#endif //NIMBLE_PENALTY_CONTACT_FORCE_CALCULATOR_H
