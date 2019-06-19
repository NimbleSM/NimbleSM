#include <iostream>
#include <vector>

int main() {
    std::cout << "Hello, world!" << '\n';
    std::vector<int> v; 
    v.push_back(10); 
    for(auto i : v) { std::cout << i << '\n'; }
}
