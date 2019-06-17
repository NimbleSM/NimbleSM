#pragma once
#include <vector>
#include "WaitAnyResult.h"
#include "mpi_err.h"

class RequestQueue
{
   public:
    std::vector<MPI_Request> pendingRequests;

   private:
    void clearAt(int index)
    {
        pendingRequests[index] = pendingRequests.back();
        pendingRequests.pop_back();
    }

   public:
    RequestQueue()                    = default;
    RequestQueue(RequestQueue const&) = delete;
    RequestQueue(RequestQueue&& queue)
      : pendingRequests(std::move(queue.pendingRequests))
    {
        queue.pendingRequests.clear();
    }

    /*
     * Adds an uninitialized request to the request queue
     * and returns a reference to that request
     */
    auto push() -> MPI_Request&
    {
        pendingRequests.emplace_back(MPI_Request{});
        return pendingRequests.back();
    }
    auto push(MPI_Request const& request) -> MPI_Request&
    {
        pendingRequests.emplace_back(request);
        return pendingRequests.back();
    }

    auto New() -> MPI_Request* { return &push(); }
    /*
     * Returns true if there are no active pendingRequests
     */
    auto empty() const -> bool { return pendingRequests.size() == 0; }
    /*
     * Returns true if there are active pendingRequests
     */
    auto has() const -> bool { return pendingRequests.size() != 0; }
    /*
     * Blocks until an active request finishes;
     * returns the request and the status and removes the request
     * from the queue of active pendingRequests
     */
    auto pop() -> MPI_Status
    {
        MPI_Status status;
        int        index;
        MPI_Waitany(
            pendingRequests.size(), pendingRequests.data(), &index, &status);
        clearAt(index);
        return status;
    }
    auto popWithSpecial(MPI_Request special, bool& popped_special) -> MPI_Status
    {
        MPI_Status status;
        int        index;
        pendingRequests.push_back(special);

        mpi_err(__FUNCTION__, ": Before Waitany", pendingRequests);
        MPI_Waitany(
            pendingRequests.size(), pendingRequests.data(), &index, &status);
        mpi_err(__FUNCTION__, ": After Waitany: ", pendingRequests);
        pendingRequests.pop_back();

        popped_special = index == pendingRequests.size();
        if (not popped_special)
        {
            clearAt(index);
        }
        return status;
    }
    auto hasMoreThan(size_t min) const -> bool
    {
        return pendingRequests.size() > min;
    }

    auto cancel_remaining() -> void
    {
        for (auto& request : pendingRequests)
        {
            MPI_Request_free(&request);
        }
        // TO-DO
    }

    auto wait_all() -> void { wait_all(MPI_STATUSES_IGNORE); }
    auto wait_all(MPI_Status* status_buffer) -> void
    {
        MPI_Waitall(pendingRequests.size(), pendingRequests.data(), status_buffer);
        pendingRequests.clear();
    }

    auto size() const noexcept -> size_t { return pendingRequests.size(); }

    ~RequestQueue() { cancel_remaining(); }
};
