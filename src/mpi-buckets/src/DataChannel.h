#pragma once
#include "RequestQueue.h"
#include "meta.h"

#include <mpi.h>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <vector>


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
    auto notify(int reciever) const -> MPI_Request
    {
        MPI_Request request;
        // MPI_Irecv(NullOf<int>(), 0, MPI_INT, reciever, tag, comm, &request);
        MPI_Isend(NullOf<int>(), 0, MPI_INT, reciever, tag, comm, &request);
        return request;
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
    auto Issend(T* dataBuffer, size_t count, int dest) const -> MPI_Request
    {
        static_assert(std::is_trivially_copyable<T>::value,
                      "Data type must be trivially copyable");
        MPI_Request request;
        MPI_Issend(
            (void*)dataBuffer, count * sizeof(T), MPI_BYTE, dest, tag, comm, &request);
        return request;
    }

    template <class T>
    auto Issend(std::vector<T> const& message, int dest) const -> MPI_Request
    {
        return Issend(message.data(), message.size(), dest);
    }

    template <class T, class T2 = T>
    auto exchange(std::vector<T> const& data, std::vector<int> const& ranks)
        -> std::vector<std::vector<T2>>
    {
        size_t                        outgoing_size = data.size();
        std::vector<size_t>           incoming_sizes(ranks.size());
        std::unordered_map<int, int>  inverseSourceRanks;
        std::vector<std::vector<int>> incomingRanks(ranks.size());
        RequestQueue                  sizeRecvQueue;
        RequestQueue                  send_queue;
        RequestQueue                  dataRecvQueue;
        inverseSourceRanks.reserve(ranks.size());

        auto dataChannel = DataChannel{comm, tag + 1};


        std::vector<std::vector<T2>> incoming(ranks.size());

        for (int i = 0; i < ranks.size(); i++)
        {
            int rank                 = ranks[i];
            inverseSourceRanks[rank] = i;
            sizeRecvQueue.push(Irecv(&incoming_sizes[i], 1, rank));
            send_queue.push(Isend(&outgoing_size, 1, rank));
            send_queue.push(dataChannel.Isend(data, rank));
        }


        while (sizeRecvQueue.has())
        {
            auto   request = sizeRecvQueue.pop();
            int    index   = inverseSourceRanks.at(request.status.MPI_SOURCE);
            size_t incoming_size = incoming_sizes[index];
            incoming[index].resize(incoming_size);
            dataRecvQueue.push(dataChannel.Irecv(
                incoming[index].data(), incoming_size, request.status.MPI_SOURCE));
        }
        dataRecvQueue.wait_all();
        send_queue.wait_all();

        return incoming;
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
