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

#ifndef NIMBLE_EXPRESSION_PARSER_H
#define NIMBLE_EXPRESSION_PARSER_H

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include <typeinfo>

#ifdef NIMBLE_HAVE_DARMA
// darma wraps headers for std containers
#include <darma.h>
#else
#include <map>
#include <string>
#endif

namespace ExpressionParsing {

struct DeletableObjectBaseClass
{  // Used by the MemoryManager class
  virtual ~DeletableObjectBaseClass(){};
};

class MemoryManager : DeletableObjectBaseClass
{
  // Keeps track of all objects created using the ``New'' method
  // When the MemoryManager is deleted, everything it manages is deleted too
 public:
  const static unsigned int BlockSize = 62;
  struct MemoryBlock
  {
    DeletableObjectBaseClass* Pointers[BlockSize];
    size_t                    count = 0;
    MemoryBlock*              next;
    MemoryBlock() : next(0) {}
    MemoryBlock(MemoryBlock* next) : next(next) {}
  };
  MemoryBlock* first;
  MemoryManager() : first(new MemoryBlock()) {}
  template <typename T, typename... In>
  T*
  New(In... inputs)
  {
    T*            pointer         = new T(inputs...);
    unsigned long insert_location = first->count;
    if (insert_location < BlockSize) {
      first->count                     = insert_location + 1;
      first->Pointers[insert_location] = pointer;
    } else {
      first              = new MemoryBlock(first);
      first->count       = 1;
      first->Pointers[0] = pointer;
    }
    return pointer;
  }
  ~MemoryManager()
  {  // Deletes everything allocated with MemoryManager::New
    MemoryBlock*               current = first;
    MemoryBlock*               next    = current->next;
    DeletableObjectBaseClass** scan    = current->Pointers;
    unsigned int               count   = current->count;
    for (unsigned int i = 0; i < count; ++i) { delete *scan++; }
    delete current;
    current = next;
    while (current) {
      scan = current->Pointers;
      next = current->next;
      for (unsigned int i = 0; i < BlockSize; ++i) { delete *scan++; }
      delete current;
      current = next;
    }
  }
};

template <typename T>
struct Evaluatable : DeletableObjectBaseClass
{  // Acts as the ``evaluation'' function (overloaded in derived classes)
  virtual operator T() const = 0;
  // Provides a way for derived classes to optimize calculation for faster
  // evaluation
  virtual Evaluatable*
  optimize(MemoryManager& m)
  {
    return this;
  }
};

template <typename T>
struct GenericConstant : Evaluatable<T>
{  // Stores a constant; when evaluated returns the constant.
  T value;
    operator T() const { return value; }
  GenericConstant() : value() {}
  GenericConstant(T value) : value(value) {}
};

template <typename T>
struct GenericReference : Evaluatable<T>
{  // Stores a pointer; when evaluated returns the value stared at the reference
  T* value;
     operator T() const { return *value; }
  GenericReference() : value(0) {}
  GenericReference(T* value) : value(value) {}
};

// This macro produces a structure which calls f(x) on the input and returns the
// result
#define MakeGenericFunction(GenericFunction, f)                                                \
  template <typename Out, typename In>                                                         \
  struct GenericFunction : Evaluatable<Out>                                                    \
  {                                                                                            \
    struct ReferenceCase : Evaluatable<Out>                                                    \
    {                                                                                          \
      In* input1;                                                                              \
          operator Out() const { return f(*input1); }                                          \
      ReferenceCase(In* in1) : input1(in1) {}                                                  \
    };                                                                                         \
    Evaluatable<In>* input1;                                                                   \
                     operator Out() const { return f(*input1); }                               \
    GenericFunction(Evaluatable<In>* in1) : input1(in1) {}                                     \
    Evaluatable<Out>*                                                                          \
    optimize(MemoryManager& m)                                                                 \
    {                                                                                          \
      Evaluatable<In>* in1 = this->input1->optimize(m);                                        \
      if (GenericConstant<In>* constant = dynamic_cast<GenericConstant<In>*>(in1)) {           \
        return m.New<GenericConstant<Out>>(f(*in1));                                           \
      } else if (GenericReference<In>* reference = dynamic_cast<GenericReference<In>*>(in1)) { \
        return m.New<ReferenceCase>(reference->value);                                         \
      }                                                                                        \
      this->input1 = in1;                                                                      \
      return this;                                                                             \
    }                                                                                          \
  };

#define MakeGenericOperator(name, operation)                                         \
  template <typename Out, typename In1, typename In2>                                \
  struct name : Evaluatable<Out>                                                     \
  {                                                                                  \
    template <typename In1_t, typename In2_t>                                        \
    struct ReferenceCase : Evaluatable<Out>                                          \
    {                                                                                \
      In1_t* input1;                                                                 \
      In2_t* input2;                                                                 \
             operator Out() const { return (*input1)operation(*input2); }            \
      ReferenceCase(In1_t* in1, In2_t* in2) : input1(in1), input2(in2) {}            \
    };                                                                               \
    template <typename In2_type>                                                     \
    struct ConstCase1 : Evaluatable<Out>                                             \
    {                                                                                \
      In1       input1;                                                              \
      In2_type* input2;                                                              \
                operator Out() const { return input1 operation(*input2); }           \
      ConstCase1(In1 in1, In2_type* in2) : input1(in1), input2(in2) {}               \
    };                                                                               \
    template <typename In1_type>                                                     \
    struct ConstCase2 : Evaluatable<Out>                                             \
    {                                                                                \
      In1_type* input1;                                                              \
      In2       input2;                                                              \
                operator Out() const { return (*input1)operation input2; }           \
      ConstCase2(In1_type* in1, In2 in2) : input1(in1), input2(in2) {}               \
    };                                                                               \
    Evaluatable<In1>* input1;                                                        \
    Evaluatable<In2>* input2;                                                        \
                      operator Out() const { return (*input1)operation(*input2); }   \
    name(Evaluatable<In1>* in1, Evaluatable<In2>* in2) : input1(in1), input2(in2) {} \
    Evaluatable<Out>*                                                                \
    optimize(MemoryManager& m)                                                       \
    {                                                                                \
      Evaluatable<In1>* input1 = this->input1->optimize(m);                          \
      Evaluatable<In2>* input2 = this->input2->optimize(m);                          \
      if (auto const1 = dynamic_cast<GenericConstant<In1>*>(input1)) {               \
        if (auto const2 = dynamic_cast<GenericConstant<In2>*>(input2))               \
          return m.New<GenericConstant<Out>>((*input1)operation(*input2));           \
        else if (auto ref2 = dynamic_cast<GenericReference<In2>*>(input2))           \
          return m.New<ConstCase1<In2>>(const1->value, ref2->value);                 \
        else                                                                         \
          return m.New<ConstCase1<Evaluatable<In2>>>(const1->value, input2);         \
      } else if (auto ref1 = dynamic_cast<GenericReference<In1>*>(input1)) {         \
        if (auto const2 = dynamic_cast<GenericConstant<In2>*>(input2))               \
          return m.New<ConstCase2<In1>>(ref1->value, const2->value);                 \
        else if (auto ref2 = dynamic_cast<GenericReference<In2>*>(input2))           \
          return m.New<ReferenceCase<In1, In2>>(ref1->value, ref2->value);           \
        else                                                                         \
          return m.New<ReferenceCase<In1, Evaluatable<In2>>>(ref1->value, input2);   \
      } else if (auto const2 = dynamic_cast<GenericConstant<In2>*>(input2))          \
        return m.New<ConstCase2<Evaluatable<In1>>>(input1, const2->value);           \
      else if (auto ref2 = dynamic_cast<GenericReference<In2>*>(input2))             \
        return m.New<ReferenceCase<Evaluatable<In1>, In2>>(input1, ref2->value);     \
      this->input1 = input1;                                                         \
      this->input2 = input2;                                                         \
      return this;                                                                   \
    }                                                                                \
  };

#define MakeGenericTwoInputFunc(name, f)                                             \
  template <typename Out, typename In1, typename In2>                                \
  struct name : Evaluatable<Out>                                                     \
  {                                                                                  \
    template <typename In1_t, typename In2_t>                                        \
    struct ReferenceCase : Evaluatable<Out>                                          \
    {                                                                                \
      In1_t* input1;                                                                 \
      In2_t* input2;                                                                 \
             operator Out() const { return f(*input1, *input2); }                    \
      ReferenceCase(In1_t* in1, In2_t* in2) : input1(in1), input2(in2) {}            \
    };                                                                               \
    template <typename In2_type>                                                     \
    struct ConstCase1 : Evaluatable<Out>                                             \
    {                                                                                \
      In1       input1;                                                              \
      In2_type* input2;                                                              \
                operator Out() const { return f(input1, *input2); }                  \
      ConstCase1(In1 in1, In2_type* in2) : input1(in1), input2(in2) {}               \
    };                                                                               \
    template <typename In1_type>                                                     \
    struct ConstCase2 : Evaluatable<Out>                                             \
    {                                                                                \
      In1_type* input1;                                                              \
      In2       input2;                                                              \
                operator Out() const { return f(*input1, input2); }                  \
      ConstCase2(In1_type* in1, In2 in2) : input1(in1), input2(in2) {}               \
    };                                                                               \
    Evaluatable<In1>* input1;                                                        \
    Evaluatable<In2>* input2;                                                        \
                      operator Out() const { return f(*input1, *input2); }           \
    name(Evaluatable<In1>* in1, Evaluatable<In2>* in2) : input1(in1), input2(in2) {} \
    Evaluatable<Out>*                                                                \
    optimize(MemoryManager& m)                                                       \
    {                                                                                \
      Evaluatable<In1>* input1 = this->input1->optimize(m);                          \
      Evaluatable<In2>* input2 = this->input2->optimize(m);                          \
      if (auto const1 = dynamic_cast<GenericConstant<In1>*>(input1)) {               \
        if (auto const2 = dynamic_cast<GenericConstant<In2>*>(input2))               \
          return m.New<GenericConstant<Out>>(f(*input1, *input2));                   \
        else if (auto ref2 = dynamic_cast<GenericReference<In2>*>(input2))           \
          return m.New<ConstCase1<In2>>(const1->value, ref2->value);                 \
        else                                                                         \
          return m.New<ConstCase1<Evaluatable<In2>>>(const1->value, input2);         \
      } else if (auto ref1 = dynamic_cast<GenericReference<In1>*>(input1)) {         \
        if (auto const2 = dynamic_cast<GenericConstant<In2>*>(input2))               \
          return m.New<ConstCase2<In1>>(ref1->value, const2->value);                 \
        else if (auto ref2 = dynamic_cast<GenericReference<In2>*>(input2))           \
          return m.New<ReferenceCase<In1, In2>>(ref1->value, ref2->value);           \
        else                                                                         \
          return m.New<ReferenceCase<In1, Evaluatable<In2>>>(ref1->value, input2);   \
      } else if (auto const2 = dynamic_cast<GenericConstant<In2>*>(input2))          \
        return m.New<ConstCase2<Evaluatable<In1>>>(input1, const2->value);           \
      else if (auto ref2 = dynamic_cast<GenericReference<In2>*>(input2))             \
        return m.New<ReferenceCase<Evaluatable<In1>, In2>>(input1, ref2->value);     \
      this->input1 = input1;                                                         \
      this->input2 = input2;                                                         \
      return this;                                                                   \
    }                                                                                \
  };

MakeGenericFunction(GenericSin, std::sin) MakeGenericFunction(GenericCos, std::cos) MakeGenericFunction(
    GenericTan,
    std::tan) MakeGenericFunction(GenericErf, std::erf) MakeGenericFunction(GenericExp, std::exp)
    MakeGenericFunction(GenericLog, std::log) MakeGenericFunction(GenericAbs, std::abs)
        MakeGenericFunction(GenericASin, std::asin) MakeGenericFunction(GenericACos, std::acos)
            MakeGenericFunction(GenericATan, std::atan) MakeGenericFunction(GenericSqrt, std::sqrt)
                MakeGenericFunction(GenericCbrt, std::cbrt) MakeGenericFunction(GenericErfc, std::erfc)
                    MakeGenericFunction(GenericCeil, std::ceil) MakeGenericFunction(GenericRound, std::round)
                        MakeGenericFunction(GenericFloor, std::floor) MakeGenericFunction(GenericLog10, std::log10)
                            MakeGenericFunction(GenericNeg, -) MakeGenericFunction(GenericNot, !)

                                MakeGenericOperator(GenericAdd, +) MakeGenericOperator(GenericMultiply, *)
                                    MakeGenericOperator(GenericDivide, /) MakeGenericOperator(GenericSubtract, -)
                                        MakeGenericOperator(GenericGreater, >) MakeGenericOperator(
                                            GenericGreaterOrEqual,
                                            >=) MakeGenericOperator(GenericLess, <)
                                            MakeGenericOperator(GenericLessOrEqual, <=)
                                                MakeGenericOperator(GenericEqual, ==)
                                                    MakeGenericOperator(GenericInequal, !=)
                                                        MakeGenericOperator(GenericAnd, &&)
                                                            MakeGenericOperator(GenericXor, xor)
                                                                MakeGenericOperator(GenericOr, ||)
                                                                    MakeGenericTwoInputFunc(GenericPow, std::pow)
                                                                        MakeGenericTwoInputFunc(GenericMod, std::fmod)

                                                                            template <typename Out>
                                                                            struct Conditional : Evaluatable<Out>
{
  Evaluatable<bool>* condition;
  Evaluatable<Out>*  IfTrue;
  Evaluatable<Out>*  IfFalse;
                     operator Out() const { return (*condition) ? (*IfTrue) : (*IfFalse); }
  Conditional(Evaluatable<bool>* condition, Evaluatable<Out>* WhenTrue, Evaluatable<Out>* WhenFalse)
      : condition(condition), IfTrue(WhenTrue), IfFalse(WhenFalse)
  {
  }
  Evaluatable<Out>*
  optimize(MemoryManager& m)
  {
    Evaluatable<bool>* condition = this->condition->optimize(m);
    Evaluatable<Out>*  IfTrue    = this->IfTrue->optimize(m);
    Evaluatable<Out>*  IfFalse   = this->IfFalse->optimize(m);
    if (auto ConstCondition = dynamic_cast<GenericConstant<bool>*>(condition)) {
      if (ConstCondition->value) {
        return IfTrue;
      } else
        return IfFalse;
    }
    this->condition = condition;
    this->IfTrue    = IfTrue;
    this->IfFalse   = IfFalse;
    return this;
  }
};

typedef Evaluatable<double>                     RealValuedExpression;
typedef GenericConstant<double>                 DoubleConstant;
typedef GenericAdd<double, double, double>      DoubleAdd;
typedef GenericSubtract<double, double, double> DoubleSubtract;
typedef GenericMultiply<double, double, double> DoubleMultiply;
typedef GenericDivide<double, double, double>   DoubleDivide;
typedef GenericPow<double, double, double>      DoublePow;
typedef GenericMod<double, double, double>      DoubleMod;
typedef GenericPow<double, double, int>         DoubleToIntegerPow;
typedef Evaluatable<bool>                       BooleanExpression;
typedef GenericConstant<bool>                   BooleanConstant;
typedef Conditional<double>                     RealValuedConditional;

struct Reader
{
  const char*  first;
  unsigned int length;
  unsigned int pos = 0;
  Reader(const char* lit, unsigned int len) : first(lit), length(len) {}
  Reader(const std::string& s) : Reader(s.c_str(), s.size()) {}
  Reader(const char* lit) : first(lit), length(std::strlen(lit)) {}
  Reader(const char* start, const char* end) : first(start), length(end - start) {}
  Reader(const Reader& r) : first(r.first), length(r.length) {}

  Reader
  sub(const char* start, const char* end)
  {
    return Reader(start, end - start);
  }

  Reader
  sub(unsigned int start, unsigned int end)
  {
    return Reader(first + start, first + end);
  }

  Reader
  sub(const char* start, unsigned int len)
  {
    return Reader(start, len);
  }

  int
  Next(char open, char close)
  {
    const char* f = this->first;
    if (f[pos] == open) {
      int depth = 1;
      while (depth != 0 && pos < length) {
        char c = f[++pos];
        if (c == open)
          depth++;
        else if (c == close)
          depth--;
      }
      return depth;
    } else
      pos++;
    return 0;
  }
  int
  Prev(char open, char close)
  {
    const char* f = this->first;
    if (f[pos] == close) {
      int depth = 1;
      while (depth != 0 && pos > 0) {
        char c = f[--pos];
        if (c == close)
          depth++;
        else if (c == open)
          depth--;
      }
      return depth;
    } else
      pos--;
    return 0;
  }

  bool
  operator>>(char c)
  {
    pos                 = 0;
    const char* f       = first;
    int         max_pos = length - 1;
    while (pos < max_pos && f[pos] != c) Next('(', ')');
    if (f[pos] == c)
      return true;
    else
      return false;
  }

  bool
  operator<<(char c)
  {
    pos           = length - 1;
    const char* f = first;
    while (pos > 0 && f[pos] != c) Prev('(', ')');
    if (f[pos] == c) return true;
    return false;
  }

  bool
  operator==(std::string& s)
  {
    return operator==(s.c_str());
  }
  bool
  operator==(const char* lit)
  {
    const char* scan = first;
    const char* last = scan + length;
    while (scan < last && *lit != '\0') {
      if (*scan == *lit) {
        scan++;
        lit++;
      } else
        return false;
    }
    return scan == last && *lit == '\0';  // Ensures that all the characters matched.
  }
  inline operator char() { return first[pos]; }
  char
  operator[](unsigned int i)
  {
    return first[i];
  }
  inline std::string
  MakeString()
  {
    return std::string(first, length);
  }
  inline bool
  IsEmpty()
  {
    return length == 0;
  }
};

bool
IsCMathFunc(Reader r);

RealValuedExpression*
GetCMathFunc(Reader r, RealValuedExpression* input, MemoryManager& m);

bool
IsDigit(char c);

bool
IsInteger(Reader r);

bool
IsNumber(Reader r);

void
ConvertStringToLowercase(std::string& s);

void
RemoveWhitespace(Reader& text);

void
RemoveRedundantParenthesis(Reader& text);

void
ApplyNecessaryFormatting(Reader& text);

class EquationContext
{
 public:
  // When the EquationContext is annihilated, the MemoryManager object will
  // annihilate everything allocated under mem.New
  MemoryManager                                mem;
  std::map<std::string, RealValuedExpression*> variables;
  void
  AddVarible(std::string name, double& ref)
  {
    variables[name] = mem.New<GenericReference<double>>(&ref);
  }
  RealValuedExpression*
  ParseEquation(std::string s)
  {
    ConvertStringToLowercase(s);
    Reader reader = Reader(s);
    return ParseEquation(reader);
  }

 private:
  RealValuedExpression*
  ParseEquation(Reader r);
  DoubleConstant*
  ParseDoubleConstant(Reader r)
  {
    if (r == "e") return mem.New<DoubleConstant>(M_E);
    if (r == "pi") return mem.New<DoubleConstant>(M_PI);
    if (r == "tau") return mem.New<DoubleConstant>(M_PI * 2);
    if (IsNumber(r)) {
      // DJL, fix for clang compiler
      // return mem.New<DoubleConstant>(strtod(r.first, NULL));
      std::string temp_string          = r.MakeString();
      double      floating_point_value = std::stod(temp_string);
      return mem.New<DoubleConstant>(floating_point_value);
    }
    return 0;
  }
  RealValuedExpression*
  ParseFunction(Reader r)
  {
    char c = *r.first;
    if (IsDigit(c) || c == '(' || c == '-' || r.length < 3) return 0;
    while (c = r.first[r.pos], c != '(' && c != ' ' && r.pos < r.length) {
      if (c != '*' && c != '/' && c != '^' && c != '+' && c != '=' && c != '<' && c != '~' && c != '>')
        r.pos++;
      else
        return 0;
    }
    if (r >> '(') {
      int p     = r.pos;
      int depth = r.Next('(', ')');
      if (depth > 0) throw std::invalid_argument("Mismatched parethesis in " + r.MakeString());
      if (r.pos < r.length - 1) return 0;
      Reader func_name = r.sub(r.first, p);
      if (IsCMathFunc(func_name)) {
        RealValuedExpression* exp = ParseEquation(r.sub(p, r.length));
        return GetCMathFunc(func_name, exp, mem);
      }
    }
    return 0;
  }
  RealValuedExpression*
  ParseOperation(Reader r)
  {
    if (r >> '?') {
      int p1    = r.pos;
      int depth = r.Next('?', ':');
      if (depth > 0) throw std::invalid_argument("Couldn't find matching : for ternary operator in " + r.MakeString());
      if (r.pos != r.length) {
        return mem.New<RealValuedConditional>(
            ParseBooleanExpression(r.sub(r.first, p1)),
            ParseEquation(r.sub(p1 + 1, r.pos)),
            ParseEquation(r.sub(r.pos + 1, r.length)));
      }
    }
    if (r << '%')
      return mem.New<DoubleMod>(ParseEquation(r.sub(r.first, r.pos)), ParseEquation(r.sub(r.pos + 1, r.length)));

    if (r << '+')
      return mem.New<DoubleAdd>(ParseEquation(r.sub(r.first, r.pos)), ParseEquation(r.sub(r.pos + 1, r.length)));

    if (r << '-') {
      // We need to ensure that we found a minus sign and not a value negation
      while (r.pos > 0 && (r.first[r.pos] == '-' || r.first[r.pos] == ' ')) r.pos--;
      if (r.pos > 0 || r.first[r.pos] != '-') {
        char c = r.first[r.pos];
        if (c != '*' && c != '/' && c != '^' && c != 'e' && c != '=' && c != '<' && c != '~' && c != '>') {
          while (r.first[r.pos] != '-') r.pos++;  // Moves it back to the last minus sign
          return mem.New<DoubleSubtract>(
              ParseEquation(r.sub(r.first, r.pos)), ParseEquation(r.sub(r.pos + 1, r.length)));
        }
      }
    }

    if (r << '*')
      return mem.New<DoubleMultiply>(ParseEquation(r.sub(r.first, r.pos)), ParseEquation(r.sub(r.pos + 1, r.length)));

    if (r << '/')
      return mem.New<DoubleDivide>(ParseEquation(r.sub(r.first, r.pos)), ParseEquation(r.sub(r.pos + 1, r.length)));

    if (r >> '^') {
      RealValuedExpression* base       = ParseEquation(r.sub(r.first, r.pos));
      Reader                power_text = r.sub(r.pos + 1, r.length);
      ApplyNecessaryFormatting(power_text);
      if (IsInteger(power_text))
        return mem.New<DoubleToIntegerPow>(base, mem.New<GenericConstant<int>>(strtol(power_text.first, NULL, 10)));
      else
        return mem.New<DoublePow>(base, ParseEquation(power_text));
    }

    if (*r.first == '-') { return mem.New<GenericNeg<double, double>>(ParseEquation(r.sub(1, r.length))); }

    return 0;
  }

  BooleanExpression*
  ParseBooleanExpression(Reader r)
  {
    ApplyNecessaryFormatting(r);
    if (r == "true") return mem.New<BooleanConstant>(true);
    if (r == "false") return mem.New<BooleanConstant>(false);

    if (r << '|')
      return mem.New<GenericOr<bool, bool, bool>>(
          ParseBooleanExpression(r.sub(r.first, r.pos)), ParseBooleanExpression(r.sub(r.pos + 1, r.length)));

    if (r << '^')
      return mem.New<GenericXor<bool, bool, bool>>(
          ParseBooleanExpression(r.sub(r.first, r.pos)), ParseBooleanExpression(r.sub(r.pos + 1, r.length)));

    if (r << '&')
      return mem.New<GenericAnd<bool, bool, bool>>(
          ParseBooleanExpression(r.sub(r.first, r.pos)), ParseBooleanExpression(r.sub(r.pos + 1, r.length)));

    if (r << '!') return mem.New<GenericNot<bool, bool>>(ParseBooleanExpression(r.sub(r.first + 1, r.length)));

    if (r << '>') {
      if (r.first[r.pos + 1] == '=')
        return mem.New<GenericGreaterOrEqual<bool, double, double>>(
            ParseEquation(r.sub(r.first, r.pos)), ParseEquation(r.sub(r.pos + 2, r.length)));
      else
        return mem.New<GenericGreater<bool, double, double>>(
            ParseEquation(r.sub(r.first, r.pos)), ParseEquation(r.sub(r.pos + 1, r.length)));
    }

    if (r << '<') {
      if (r.first[r.pos + 1] == '=')
        return mem.New<GenericLessOrEqual<bool, double, double>>(
            ParseEquation(r.sub(r.first, r.pos)), ParseEquation(r.sub(r.pos + 2, r.length)));
      else
        return mem.New<GenericLess<bool, double, double>>(
            ParseEquation(r.sub(r.first, r.pos)), ParseEquation(r.sub(r.pos + 1, r.length)));
    }

    /*if(r << '~') {
            if(r.first[r.pos + 1] == '=') return
    mem.New<DoubleComparison>(ApproxEqual, ParseEquation(r.sub(r.pos + 2,
    r.length)), ParseEquation(r.sub(r.first, r.pos))); else return
    mem.New<DoubleComparison>(ApproxEqual, ParseEquation(r.sub(r.pos + 1,
    r.length)), ParseEquation(r.sub(r.first, r.pos)));
    }*/

    if (r << '=') {
      if (r.first[r.pos - 1] == '=')
        return mem.New<GenericEqual<bool, double, double>>(
            ParseEquation(r.sub(r.pos + 1, r.length)), ParseEquation(r.sub(r.first, r.pos - 1)));
    }
    throw std::invalid_argument("Unable to parse \"" + r.MakeString() + "\" as boolean");
  }
};

// To make a copy of a BoundaryConditionFunctor, create a new one directly from
// the string that was originally used to parse it (i.e,
// BoundaryConditionFunctor.equation)
struct BoundaryConditionFunctor
{
  double           x = 0, y = 0, z = 0, t = 0;
  EquationContext* context;
  // The Computable Expression Object
  RealValuedExpression* expression = nullptr;
  std::string           equation;
  BoundaryConditionFunctor() : context(0) {}  // Sets context to nil
  BoundaryConditionFunctor(const std::string& equation)
  {
    if (equation == "") {
      // Acts as default constructor if given an empty string
      context = 0;
      return;
    }
    try {
      context = new EquationContext();
      // Usage: context->AddVarible(name, reference);
      context->AddVarible("x", x);  // Passes a reference to x
      context->AddVarible("y", y);  // Passes a reference to y
      context->AddVarible("z", z);  // Passes a reference to z
      context->AddVarible("t", t);  // Passes a reference to t
      // Creates the Computable Expression Object
      expression = context->ParseEquation(equation);
      // Optimizes the Computable Expression Object
      expression     = expression->optimize(context->mem);
      this->equation = equation;  // Stores a copy of the equation
    } catch (std::invalid_argument arg) {
      if (context) delete context;
      context = 0;
      throw arg;
    }
  }
  BoundaryConditionFunctor(const BoundaryConditionFunctor& other) : BoundaryConditionFunctor(other.equation) {}
  BoundaryConditionFunctor&
  operator=(const BoundaryConditionFunctor& other)
  {
    return *this = other.equation;
  }
  BoundaryConditionFunctor&
  operator=(const std::string& equation)
  {
    EquationContext* old_context = context;
    try {
      context = new EquationContext();
      context->AddVarible("x", x);
      context->AddVarible("y", y);
      context->AddVarible("z", z);
      context->AddVarible("t", t);
      expression     = context->ParseEquation(equation);
      expression     = expression->optimize(context->mem);
      this->equation = equation;
      if (old_context) delete old_context;
      return *this;
    } catch (std::invalid_argument arg) {
      if (context) delete context;
      context = old_context;
      throw arg;
    }
  }
  inline double
  eval()
  {
    if (!expression) {
      throw std::invalid_argument(
          "Error in BoundaryConditionFunctor::eva(), expression pointer is "
          "null.");
    }
    return *expression;
  }
  inline double
  eval(double x, double y, double z, double t)
  {
    this->x = x;
    this->y = y;
    this->z = z;
    this->t = t;
    if (!expression) {
      throw std::invalid_argument(
          "Error in BoundaryConditionFunctor::eva(), expression pointer is "
          "null.");
    }
    return *expression;
  }
  inline double
  operator()(double x, double y, double z, double t)
  {
    this->x = x;
    this->y = y;
    this->z = z;
    this->t = t;
    if (!expression) {
      throw std::invalid_argument(
          "Error in BoundaryConditionFunctor::eva(), expression pointer is "
          "null.");
    }
    return *expression;  // Evaluates the Computable Expression Object
  }
  ~BoundaryConditionFunctor()
  {
    if (context) delete context;
    // RealValuedExpression* expression is managed by the EquationContext,
    // And will be destroyed when "context" is deleted.
  }
};
}  // namespace ExpressionParsing

#endif
