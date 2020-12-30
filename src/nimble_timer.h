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

#ifndef NIMBLE_TIMER_H
#define NIMBLE_TIMER_H

#include <chrono>

#ifdef NIMBLE_HAVE_DARMA
  #include <darma.h>
#else
  #include <string>
  #include <map>
#endif

namespace nimble {

  class TimeKeeper {

  public:

    TimeKeeper() : total_time_(0.0) {}

    template<typename ArchiveType>
    void serialize(ArchiveType& ar) {
      ar | total_time_;
    }

    inline void Start() {
      start_time_ = std::chrono::high_resolution_clock::now();
    }

    inline void Stop() {
      using std::chrono::duration;
      using std::chrono::duration_cast;
      end_time_ = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> time_increment(0.0);
      if (end_time_ > start_time_) {
        time_increment = duration_cast<duration<double>>(end_time_ - start_time_);
      }
      total_time_ += time_increment.count();
    }

    inline double GetElapsedTime() const {
      return total_time_;
    }

  private:

    std::chrono::time_point<std::chrono::high_resolution_clock> start_time_;
    std::chrono::time_point<std::chrono::high_resolution_clock> end_time_;
    double total_time_;
  };

  class Timer {

  public:

    Timer()= default;

    ~Timer()= default;

    template<typename ArchiveType>
    void serialize(ArchiveType& ar) {
      ar | timers_;
    }

    inline void Start(const std::string& name) {
      TimeKeeper& time_keeper = timers_[name];
      time_keeper.Start();
    }

    inline void Stop(const std::string& name) {
      TimeKeeper& time_keeper = timers_[name];
      time_keeper.Stop();
    }

    inline double ElapsedTime(const std::string& name) const {
      return timers_.at(name).GetElapsedTime();
    }

    std::map<std::string, TimeKeeper> timers_;
  };
} // namespace nimble

#endif // NIMBLE_TIMER_H
