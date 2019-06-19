#pragma once
#include <vector>

struct KokkosMock
{
    std::vector<double> values;
    auto                extent(int) const -> size_t { return values.size(); }
    auto operator()(int index) const -> double const& { return values[index]; }
};
