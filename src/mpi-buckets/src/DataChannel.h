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
    auto Iawait(int sender) const -> MPI_Request {
        return Irecv(NullOf<int>(), 0, sender); 
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
};
