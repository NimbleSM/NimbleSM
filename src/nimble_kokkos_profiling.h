#pragma once

#include <nimble_kokkos_defs.h>

#include <string>

namespace nimble_kokkos {

class ProfilingTimer
{
 public:
  inline ProfilingTimer() = default;

  inline void
  push_region(const std::string& profiling_region_name)
  {
    Kokkos::Profiling::pushRegion("NimbleSM: " + profiling_region_name);
    timer.reset();
  }

  inline double
  pop_region_and_report_time() const
  {
    double t = timer.seconds();
    Kokkos::Profiling::popRegion();
    return t;
  }

 private:
  Kokkos::Timer timer;
};

}  // namespace nimble_kokkos