#include "parser.hpp"

#include <cstdio>
#include <iostream>

template <class... T>
void f()
{
    puts(__PRETTY_FUNCTION__);
}
struct CF
{
    template <class... T>
    auto operator()(T&&... t) -> void
    {
        f<decltype(t)...>();
    }
} cf;
void test_typespec()
{
    variant<int, double> a{10};
    visit(cf, a);
    visit(cf, (const variant<int, double>&)a);
    visit(cf, std::move(a));
}
int main() { test_typespec(); }
