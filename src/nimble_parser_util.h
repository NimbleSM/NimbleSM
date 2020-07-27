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

#ifndef SRC_NIMBLE_PARSER_UTIL_H_
#define SRC_NIMBLE_PARSER_UTIL_H_

#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace nimble {

namespace {

inline std::vector<std::string> tokenize_string_with_separator(const std::string &s, const char sep) {
  std::vector<std::string> tokens;
  std::stringstream ss(s);
  std::string token;
  while (std::getline(ss >> std::ws, token, sep)) {
    tokens.push_back(token);
  }
  return tokens;
}

}

inline std::vector<std::string> tokenize_string(const std::string& s) {
  std::vector<std::string> tokens;
  auto quote_delimited_segments = tokenize_string_with_separator(s, '"');
  auto is_quoted_segment = [](int i) -> bool {
    return (i % 2) == 1;
  };
  auto quote = quote_delimited_segments.begin();
  auto quote_end = quote_delimited_segments.end();
  int i = 0;
  for (; quote != quote_end; ++quote, ++i) {
    if (is_quoted_segment(i)) {
      tokens.push_back(*quote);
    } else {
      auto segment_tokens = tokenize_string_with_separator(*quote, ' ');
      tokens.insert(tokens.end(), segment_tokens.begin(), segment_tokens.end());
    }
  }
  return tokens;
}

inline double string_to_double(const std::string& s) {
  std::stringstream ss;
  ss << s;
  double val = std::numeric_limits<double>::max();
  ss >> val;
  return val;
}

}

#endif /* EXTRAS_SRC_NIMBLE_PARSER_UTIL_H_ */
