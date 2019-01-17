#pragma once
#include <vector>
#include "WaitAnyResult.h"

class RequestQueue
{
    std::vector<MPI_Request> requests;

    auto clearAt(int index) -> MPI_Request
    {
        MPI_Request request = requests[index];
        requests[index]     = requests.back();
        requests.pop_back();
        return request;
    }

   public:
    /*
     * Adds an uninitialized request to the request queue
     * and returns a reference to that request
     */
    auto push() -> MPI_Request&
    {
        requests.emplace_back(MPI_Request{});
        return requests.back();
    }
    auto push(MPI_Request const& request) -> MPI_Request&
    {
        requests.emplace_back(request);
        return requests.back();
    }

    auto New() -> MPI_Request* { return &push(); }
    /*
     * Returns true if there are no active requests
     */
    auto empty() const -> bool { return requests.size() == 0; }
    /*
     * Returns true if there are active requests
     */
    auto has() const -> bool { return requests.size() != 0; }
    /*
     * Blocks until an active request finishes;
     * returns the request and the status and removes the request
     * from the queue of active requests
     */
    auto pop() -> WaitAnyResult
    {
        auto result = WaitAnyResult{};
        int  index;
        MPI_Waitany(requests.size(), requests.data(), &index, &result.status);
        result.request = clearAt(index);
        return result;
    }
};
