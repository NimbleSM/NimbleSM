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

#ifndef NIMBLE_MACROS_H
#define NIMBLE_MACROS_H

#include <iostream>  // For std::cerr
#include <pmmintrin.h> // SSE intrinsics
#include <xmmintrin.h> // SSE intrinsics

#define NIMBLE_ASSERT_IMPL(cond, msg, ...)                              \
  do {                                                                  \
    if (!(cond)) {                                                      \
      std::cerr << #cond " NIMBLE_ASSERT failed at ";                   \
      std::cerr << __FILE__ << " +" << __LINE__ << "\n" << msg << '\n'; \
      abort();                                                          \
    }                                                                   \
  } while (0)

#define NIMBLE_PANIC_IMPL(cond, msg, ...)                               \
  do {                                                                  \
    if ((cond)) {                                                       \
      std::cerr << #cond " NIMBLE_PANIC condition at ";                 \
      std::cerr << __FILE__ << " +" << __LINE__ << "\n" << msg << '\n'; \
      abort();                                                          \
    }                                                                   \
  } while (0)

#define NIMBLE_ABORT_IMPL(msg, ...)                                   \
  do {                                                                \
    std::cerr << " NIMBLE_ABORT statement at ";                       \
    std::cerr << __FILE__ << " +" << __LINE__ << "\n" << msg << '\n'; \
    abort();                                                          \
  } while (0)

#define NIMBLE_TRACE_IMPL(msg, ...)                                   \
  do {                                                                \
    std::cerr << "********** NIMBLE_TRACE at ";                       \
    std::cerr << __FILE__ << " +" << __LINE__ << "\n" << msg << '\n'; \
  } while (0)

#define NIMBLE_DUMP_IMPL(msg, ...) \
  do {                             \
    std::cerr << msg;              \
  } while (0)

#define NIMBLE_TRAP_FPE_IMPL(...)                                                                                      \
  do {                                                                                                              \
    _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);                                                                     \
    _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);                                                             \
    _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~(_MM_MASK_INVALID | _MM_MASK_DIV_ZERO | _MM_MASK_OVERFLOW)); \
  } while (0)

#define NIMBLE_ASSERT(...) NIMBLE_ASSERT_IMPL(__VA_ARGS__, "")
#define NIMBLE_PANIC(...) NIMBLE_PANIC_IMPL(__VA_ARGS__, "")
#define NIMBLE_ABORT(...) NIMBLE_ABORT_IMPL(__VA_ARGS__, "")
#define NIMBLE_TRACE(...) NIMBLE_TRACE_IMPL(__VA_ARGS__, "")
#define NIMBLE_DUMP(...) NIMBLE_DUMP_IMPL(__VA_ARGS__, "")

#ifndef NIMBLE_DEBUG
#define NIMBLE_EXPECT(...)
#else
#define NIMBLE_EXPECT(...) NIMBLE_ASSERT(__VA_ARGS__)
#endif

#define NIMBLE_ALWAYS_ASSERT(cond) NIMBLE_ASSERT(cond)
#define NIMBLE_ALWAYS_ASSERT_VERBOSE(cond, msg) NIMBLE_ASSERT(cond, msg)
#define NIMBLE_DEBUG_ASSERT(cond) NIMBLE_EXPECT(cond)
#define NIMBLE_DEBUG_ASSERT_VERBOSE(cond, msg) NIMBLE_EXPECT(cond, msg)
#define NIMBLE_TRAP_FPE(...) NIMBLE_TRAP_FPE_IMPL(__VA_ARGS__)

#endif  // NIMBLE_MACROS_H
