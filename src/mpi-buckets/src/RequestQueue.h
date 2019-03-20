#pragma once
#include <vector>
#include "WaitAnyResult.h"

class RequestQueue
{
    std::vector<MPI_Request> pendingRequests;

    auto clearAt(int index) -> MPI_Request
    {
        MPI_Request request    = pendingRequests[index];
        pendingRequests[index] = pendingRequests.back();
        pendingRequests.pop_back();
        return request;
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
    auto pop() -> WaitAnyResult
    {
        auto result = WaitAnyResult{};
        int  index;
        MPI_Waitany(
            pendingRequests.size(), pendingRequests.data(), &index, &result.status);
        result.request = clearAt(index);
        return result;
    }
    auto popWithIndex() -> std::pair<int, WaitAnyResult>
    {
        auto result = WaitAnyResult{};
        int  index;
        MPI_Waitany(
            pendingRequests.size(), pendingRequests.data(), &index, &result.status);
        result.request = clearAt(index);
        return {index, result};
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
