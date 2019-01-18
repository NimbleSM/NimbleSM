#pragma once
#include <mpi.h>
#include <type_traits>
#include <vector>
#include "meta.h"

struct DataChannel
{
    constexpr static auto any_source = MPI_ANY_SOURCE;

    MPI_Comm comm;
    int      tag;

    auto commRank() const -> int
    {
        int rank;
        MPI_Comm_rank(comm, &rank);
        return rank;
    }
    auto commSize() const -> int
    {
        int size;
        MPI_Comm_size(comm, &size);
        return size;
    }
    auto Iawait(int sender) const -> MPI_Request
    {
        return Irecv(NullOf<int>(), 0, sender);
    }
    /**
     * @brief Fulfill a pending await on another rank
     *
     * @param reciever
     * @return void
     */
    auto notify(int reciever) const -> void
    {
        MPI_Send(NullOf<int>(), 0, MPI_INT, reciever, tag, comm);
    }
    template <class T>
    auto Irecv(T* dataBuffer, size_t count, int source) const -> MPI_Request
    {
        static_assert(std::is_trivially_copyable<T>::value,
                      "Data type must be trivially copyable");
        MPI_Request request;
        MPI_Irecv(
            (void*)dataBuffer, count * sizeof(T), MPI_BYTE, source, tag, comm, &request);
        return request;
    }
    template <class T>
    auto Irecv(std::vector<T>& message, int source) const -> MPI_Request
    {
        return Irecv(message.data(), message.size(), source);
    }
    template <class T>
    auto recv(T* dataBuffer, size_t count, int source) const -> void
    {
        MPI_Recv(
            dataBuffer, count * sizeof(T), MPI_BYTE, source, tag, comm, nullptr);
    }

    template <class T>
    auto Isend(T* dataBuffer, size_t count, int dest) const -> MPI_Request
    {
        static_assert(std::is_trivially_copyable<T>::value,
                      "Data type must be trivially copyable");
        MPI_Request request;
        MPI_Isend(
            (void*)dataBuffer, count * sizeof(T), MPI_BYTE, dest, tag, comm, &request);
        return request;
    }

    template <class T>
    auto Isend(std::vector<T> const& message, int dest) const -> MPI_Request
    {
        return Isend(message.data(), message.size(), dest);
    }

    template <class T>
    auto IrecvAny(T* dataBuffer, size_t count) const -> MPI_Request
    {
        static_assert(std::is_trivially_copyable<T>::value,
                      "Data type must be trivially copyable");
        MPI_Request request;
        MPI_Irecv((void*)dataBuffer,
                  count * sizeof(T),
                  MPI_BYTE,
                  any_source,
                  tag,
                  comm,
                  &request);
        return request;
    }
    template <class T>
    auto IrecvAny(std::vector<T>& message) const -> MPI_Request
    {
        return IrecvAny(message.data(), message.size());
    }
    template <class F, class... ExtraArgs>
    auto onChildRanks(F&& func, ExtraArgs&&... args) const -> void
    {
        int const       num_ranks = commSize();
        long long const my_rank   = commRank();
        auto const      child_0   = (my_rank + 1) * 2 - 1;
        auto const      child_1   = (my_rank + 1) * 2;
        if (child_0 < num_ranks)
        {
            func(args..., *this, (int)child_0);
        }
        if (child_1 < num_ranks)
        {
            func(args..., *this, (int)child_1);
        }
    }
    template <class F, class... ExtraArgs>
    auto onParentRanks(F&& func, ExtraArgs&&... args) const -> void
    {
        int const my_rank = commRank();
        if (my_rank != 0)
        {
            const int parent = (my_rank - 1) / 2;
            func(args..., *this, parent);
        }
    }
    auto isFromParent(MPI_Status const& status) const -> bool
    {
        return isFromParent(status.MPI_SOURCE);
    }
    auto isFromChild(MPI_Status const& status) const -> bool
    {
        return isFromChild(status.MPI_SOURCE);
    }
    auto isFromParent(int sender) const -> bool
    {
        int const parent = (commRank() - 1) / 2;
        return sender == parent;
    }
    auto isFromChild(int sender) const -> bool
    {
        return (sender - 1) / 2 == commRank();
    }
    auto isRoot() const -> bool { return commRank() == 0; }
    auto countChildren() const -> int
    {
        int count = 0;
        onChildRanks([&](DataChannel const&, int) { ++count; });
        return count;
    }
};
