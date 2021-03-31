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

#include "nimble_expression_parser.h"

bool
ExpressionParsing::IsCMathFunc(ExpressionParsing::Reader r)
{
  if (r == "sin" || r == "cos" || r == "tan" || r == "erf" || r == "exp" ||
      r == "log" || r == "abs" || r == "asin" || r == "acos" || r == "atan" ||
      r == "sqrt" || r == "cbrt" || r == "erfc" || r == "ceil" ||
      r == "round" || r == "floor" || r == "log10")
    return true;
  return false;
}

ExpressionParsing::RealValuedExpression*
ExpressionParsing::GetCMathFunc(
    ExpressionParsing::Reader                r,
    ExpressionParsing::RealValuedExpression* input,
    ExpressionParsing::MemoryManager&        m)
{
  if (r == "sin") return m.New<GenericSin<double, double>>(input);
  if (r == "cos") return m.New<GenericCos<double, double>>(input);
  if (r == "tan") return m.New<GenericTan<double, double>>(input);
  if (r == "erf") return m.New<GenericErf<double, double>>(input);
  if (r == "exp") return m.New<GenericExp<double, double>>(input);
  if (r == "log") return m.New<GenericLog<double, double>>(input);
  if (r == "abs") return m.New<GenericAbs<double, double>>(input);
  if (r == "asin") return m.New<GenericASin<double, double>>(input);
  if (r == "acos") return m.New<GenericACos<double, double>>(input);
  if (r == "atan") return m.New<GenericATan<double, double>>(input);
  if (r == "sqrt") return m.New<GenericSqrt<double, double>>(input);
  if (r == "cbrt") return m.New<GenericCbrt<double, double>>(input);
  if (r == "erfc") return m.New<GenericErfc<double, double>>(input);
  if (r == "ceil") return m.New<GenericCeil<double, double>>(input);
  if (r == "round") return m.New<GenericRound<double, double>>(input);
  if (r == "floor") return m.New<GenericFloor<double, double>>(input);
  if (r == "log10") return m.New<GenericLog10<double, double>>(input);
  return 0;
}

bool
ExpressionParsing::IsDigit(char c)
{
  return c >= '0' && c <= '9';
}

bool
ExpressionParsing::IsInteger(ExpressionParsing::Reader r)
{
  if (r.IsEmpty()) return false;
  int len = r.length;
  if (len == 1 && (r[0] == '-' || r[0] == '+')) return false;
  if (r[0] != '-' && r[0] != '+' && !IsDigit(r[0])) return false;
  for (int i = 1; i < len; i++)
    if (!IsDigit(r[i])) return false;
  return true;
}

bool
ExpressionParsing::IsNumber(ExpressionParsing::Reader r)
{
  if (r.length == 0) return false;
  r.pos        = 0;
  bool decimal = false;
  bool digits  = false;
  if (r == '-' || r == '+') r.pos++;
  while (r.pos < r.length) {
    char c = r;
    r.pos++;
    if (IsDigit(c)) {
      digits = true;
      continue;
    }
    if (c == '.' && !decimal) {
      decimal = true;
      continue;
    }
    if (c == 'e' && digits) { return IsInteger(r.sub(r.pos, r.length)); }
    return false;
  }
  return true;
}

void
ExpressionParsing::ConvertStringToLowercase(std::string& s)
{
  int   len = s.size();
  char* c   = (char*)s.c_str();
  for (int i = 0; i < len; i++) {
    if (c[i] >= 'A' && c[i] <= 'Z') c[i] += 32;
  }
}

void
ExpressionParsing::RemoveWhitespace(ExpressionParsing::Reader& text)
{
  const char* first = text.first;
  const char* last  = first + text.length - 1;
  if (*first != ' ' && *last != ' ') return;
  while (*first == ' ' && first < last) first++;
  while (*last == ' ' && last > first) last--;
  text.first  = first;
  text.length = (last - first + 1);
  return;
}

void
ExpressionParsing::RemoveRedundantParenthesis(ExpressionParsing::Reader& text)
{
  text.pos = 0;
  while (text == '(') {
    text.Next('(', ')');  // Jumps to the closing parenthesis
    if (text.pos == text.length - 1) {
      text.first++;
      text.length -= 2;
      text.pos = 0;
    } else {
      text.pos = 0;
      break;
    }
  }
}

void
ExpressionParsing::ApplyNecessaryFormatting(ExpressionParsing::Reader& text)
{
  if (text.IsEmpty()) return;
  RemoveWhitespace(text);
  RemoveRedundantParenthesis(text);
}

ExpressionParsing::RealValuedExpression*
ExpressionParsing::EquationContext::ParseEquation(ExpressionParsing::Reader r)
{
  ApplyNecessaryFormatting(r);
  if (r.IsEmpty()) throw std::invalid_argument("Can't parse empty string");

  RealValuedExpression* output = ParseDoubleConstant(r);
  if (output) return output;

  output = ParseOperation(r);
  if (output) return output;

  output = ParseFunction(r);
  if (output) return output;

  auto it = variables.find(r.MakeString());
  if (it != variables.end()) return it->second;
  throw std::invalid_argument("Unable to parse \"" + r.MakeString() + "\"");
}
