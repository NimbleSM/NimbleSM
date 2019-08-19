#include "mpi/DataChannel.hpp"
#include <array>

template<size_t n> 
struct MultiChannel {
    MPI_Comm comm; 
    std::array<int, n> tags; 
    constexpr auto get(int i) const noexcept -> DataChannel {
        return DataChannel{comm, tags[i]}; 
    }
    template<size_t i> 
    constexpr auto get() const noexcept -> DataChannel {
        static_assert(i < n, "The index must be smaller than the number of channels"); 
        return DataChannel{comm, tags[i]};
    }
    template<size_t... i> 
    constexpr auto sub() const noexcept -> MultiChannel<sizeof...(i)> {
        return { comm, tags[i]...}; 
    };
};
