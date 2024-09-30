#ifndef NIMBLE_CONTACT_FORCE_CALCULATOR_H
#define NIMBLE_CONTACT_FORCE_CALCULATOR_H

#include <array>

namespace nimble {

  class ContactForceCalculator {
  public:
    virtual ~ContactForceCalculator() = default;

    virtual void ComputeContactForce(double penalty, double gap, const double normal[3], double mass1, std::array<double, 3> a1, std::array<double, 3> a2, std::array<double, 3> N2, std::array<double, 3>& contact_force) const = 0;
  };

} // namespace nimble

#endif //NIMBLE_CONTACT_FORCE_CALCULATOR_H
