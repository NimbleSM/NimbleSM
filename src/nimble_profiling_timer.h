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

#ifndef NIMBLESM_NIMBLE_PROFILING_TIMER_H
#define NIMBLESM_NIMBLE_PROFILING_TIMER_H

#include <memory>
#include <stack>
#include <string>

#include "nimble_timer.h"

#ifdef NIMBLE_HAVE_KOKKOS
#include "Kokkos_Timer.hpp"
#endif

namespace nimble {

class ProfilingTimer
{
 public:
  inline ProfilingTimer()
  {
#ifndef NIMBLE_HAVE_KOKKOS
    n_timer_ = std::unique_ptr<nimble::Timer>(new nimble::Timer());
#endif
  }

  ~ProfilingTimer() = default;

  inline void
  push_region(const std::string& profiling_region_name)
  {
    current_region_.push(profiling_region_name);
#ifdef NIMBLE_HAVE_KOKKOS
    Kokkos::Profiling::pushRegion("NimbleSM: " + profiling_region_name);
    k_timer_.reset();
#else
    n_timer_->Start(profiling_region_name);
#endif
  }

  inline double
  pop_region_and_report_time()
  {
    double t = 0.0;
#ifdef NIMBLE_HAVE_KOKKOS
    t = k_timer_.seconds();
    Kokkos::Profiling::popRegion();
#else
    t = n_timer_->StopReport(current_region_.top());
#endif
    current_region_.pop();
    return t;
  }

 private:
  std::stack<std::string> current_region_;

  std::unique_ptr<nimble::Timer> n_timer_ = nullptr;
#ifdef NIMBLE_HAVE_KOKKOS
  Kokkos::Timer k_timer_;
#endif
};

}  // namespace nimble

#endif  // NIMBLESM_NIMBLE_PROFILING_TIMER_H
