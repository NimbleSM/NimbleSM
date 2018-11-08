/*
//@HEADER
// ************************************************************************
//
//                NimbleSM
//               Copyright 2018
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

#ifndef NIMBLE_LIGHTWEIGHT_STOPWATCH_H
#define NIMBLE_LIGHTWEIGHT_STOPWATCH_H

#include <chrono>
#include <ctime>
#include <iomanip>
#include <sstream>
#include <thread>
#include <string>

namespace nimble
{
namespace quanta
{
struct stopwatch
{
  typedef decltype(std::chrono::high_resolution_clock::now()) time;
  typedef std::chrono::duration<double, std::ratio<1, 1>> seconds;
  typedef std::chrono::duration<double, std::milli> milliseconds;
  typedef std::chrono::duration<double, std::micro> microseconds;
  typedef decltype((time{} - time{}).count()) nanoseconds;
  typedef decltype((time{} - time{})) duration;
  time _start = std::chrono::high_resolution_clock::now();
  nanoseconds age_nano() { return (std::chrono::high_resolution_clock::now() - _start).count(); }
  double age_micro()
  {
    return microseconds(std::chrono::high_resolution_clock::now() - _start).count();
  }
  double age_milli()
  {
    return milliseconds(std::chrono::high_resolution_clock::now() - _start).count();
  }
  double age() { return seconds(std::chrono::high_resolution_clock::now() - _start).count(); }
  void reset() { _start = std::chrono::high_resolution_clock::now(); }
  template<class Function, class... Args>
  static double timefunction(Function&& f, Args&&... inputs)
  {
    stopwatch s;
    s.reset();
    f(std::forward<Args>(inputs)...);
    return s.age();
  }
  static std::string get_nanosecond_timestamp()
  {
    auto epoch_duration      = std::chrono::high_resolution_clock::now().time_since_epoch();
    auto duration_in_seconds = std::chrono::duration_cast<std::chrono::seconds>(epoch_duration);
    auto nanoseconds =
        (std::chrono::duration_cast<std::chrono::nanoseconds>(epoch_duration - duration_in_seconds)
             .count());
    std::time_t t = duration_in_seconds.count();
    char buffer[100]{};
    size_t end;
    if ((end = std::strftime(buffer, 100, "%Y.%m.%e.%H.%M.%S.", std::localtime(&t))))
    {
      buffer[end]     = '0' + (nanoseconds / 100000000 % 10);
      buffer[end + 1] = '0' + (nanoseconds / 10000000 % 10);
      buffer[end + 2] = '0' + (nanoseconds / 1000000 % 10);
      buffer[end + 3] = '0' + (nanoseconds / 100000 % 10);
      buffer[end + 4] = '0' + (nanoseconds / 10000 % 10);
      buffer[end + 5] = '0' + (nanoseconds / 1000 % 10);
      buffer[end + 6] = '0' + (nanoseconds / 100 % 10);
      buffer[end + 7] = '0' + (nanoseconds / 10 % 10);
      buffer[end + 8] = '0' + (nanoseconds % 10);
      buffer[end + 9] = '\0';
    }
    return buffer;
  }
  static std::string get_microsecond_timestamp()
  {
    auto epoch_duration      = std::chrono::high_resolution_clock::now().time_since_epoch();
    auto duration_in_seconds = std::chrono::duration_cast<std::chrono::seconds>(epoch_duration);
    auto microseconds =
        (std::chrono::duration_cast<std::chrono::microseconds>(epoch_duration - duration_in_seconds)
             .count());
    std::time_t t = duration_in_seconds.count();
    char buffer[100]{};
    size_t end;
    if ((end = std::strftime(buffer, 100, "%Y.%m.%d.%H.%M.%S.", std::localtime(&t))))
    {
      buffer[end]     = '0' + (microseconds / 100000 % 10);
      buffer[end + 1] = '0' + (microseconds / 10000 % 10);
      buffer[end + 2] = '0' + (microseconds / 1000 % 10);
      buffer[end + 3] = '0' + (microseconds / 100 % 10);
      buffer[end + 4] = '0' + (microseconds / 10 % 10);
      buffer[end + 5] = '0' + (microseconds % 10);
      buffer[end + 6] = '\0';
    }
    return buffer;
  }
  static std::string get_millisecond_timestamp()
  {
    auto epoch_duration      = std::chrono::high_resolution_clock::now().time_since_epoch();
    auto duration_in_seconds = std::chrono::duration_cast<std::chrono::seconds>(epoch_duration);
    auto milliseconds =
        (std::chrono::duration_cast<std::chrono::milliseconds>(epoch_duration - duration_in_seconds)
             .count());
    std::time_t t = duration_in_seconds.count();
    char buffer[100]{};
    size_t end;
    if ((end = std::strftime(buffer, 100, "%Y.%m.%e.%H.%M.%S.", std::localtime(&t))))
    {
      size_t milliseconds = (std::chrono::duration_cast<std::chrono::milliseconds>(
                                 std::chrono::high_resolution_clock::now().time_since_epoch())
                                 .count());
      buffer[end]         = '0' + (milliseconds / 100 % 10);
      buffer[end + 1]     = '0' + (milliseconds / 10 % 10);
      buffer[end + 2]     = '0' + (milliseconds % 10);
      buffer[end + 3]     = '\0';
    }
    return buffer;
  }
};
}   // namespace quanta
}   // namespace nimble

#endif // NIMBLE_LIGHTWEIGHT_STOPWATCH_H
