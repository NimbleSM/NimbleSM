#include <iostream>
#include <vector>
#include <algorithm>
#include <utility>
#include <numeric>
#include "View.h"
int main() {
    std::vector<int> v(100);
    auto view = ViewDeets::of(v); 
    static_assert(noexcept(view.begin()),    "View.begin() expected to be noexcept");
    static_assert(noexcept(view.end()),      "View.end() expected to be noexcept"); 
    static_assert(noexcept(view.advance()),  "View.advance() expected to be noexcept");
    static_assert(noexcept(view.head()),     "View.head() expected to be noexcept");
    static_assert(noexcept(view.tail()),     "View.tail() expected to be noexcept");


    static_assert(noexcept(view.size()), "View.size() expected to be noexcept");
    std::iota(view.begin(), view.end(), 0); 
    for(int i : view) {
        std::cout << i << std::endl; 
    }
    return 0;
}
