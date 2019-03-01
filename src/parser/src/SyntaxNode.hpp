#pragma once
#include <memory>
#include <string>
#include <vector>
#include "util/variant.hpp"

class BinaryOperator;

template <template <class...> class V>
class SyntaxTypes
{
   public:
    using type = V<std::string, BinaryOperator>;
};
template <template <class...> class V>
using STypes = typename SyntaxTypes<V>::type;

class SyntaxNode;
class BinaryOperator
{
    enum Operation { Add, Sub, Mul, Div };

   public:
    Operation                   op;
    std::unique_ptr<SyntaxNode> left;
    std::unique_ptr<SyntaxNode> right;
};

class SyntaxNode
{
    variant<std::string, BinaryOperator> node;
};
