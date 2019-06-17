#pragma once
#include <tuple>
#include <vector>
#include "WaitAnyResult.hpp"
#include "mpi_err.hpp"

template <class Tag>
class TaggedRequestQueue
{
   public:
    std::vector<MPI_Request> pendingRequests;
    std::vector<Tag>         tags;

   private:
    auto popAt(int index) -> Tag
    {
        pendingRequests[index] = pendingRequests.back();
        Tag tag                = tags[index];
        tags[index]            = tags.back();
        pendingRequests.pop_back();
        tags.pop_back();
        return tag;
    }

   public:
    TaggedRequestQueue()                           = default;
    TaggedRequestQueue(TaggedRequestQueue const&)  = delete;
    TaggedRequestQueue(TaggedRequestQueue&& queue) = default;

    auto push(MPI_Request const& request, Tag const& tag) -> void
    {
        pendingRequests.emplace_back(request);
        tags.emplace_back(tag);
    }
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
    auto pop() -> std::pair<MPI_Status, Tag>
    {
        MPI_Status status;
        int        index;
        MPI_Waitany(
            pendingRequests.size(), pendingRequests.data(), &index, &status);
        return {status, popAt(index)};
    }
    auto popWithSpecial(MPI_Request special, Tag special_tag, bool& popped_special)
        -> std::pair<MPI_Status, Tag>
    {
        MPI_Status status;
        int        index;
        pendingRequests.push_back(special);

        mpi_err(__FUNCTION__, ": Before Waitany", pendingRequests);
        MPI_Waitany(
            pendingRequests.size(), pendingRequests.data(), &index, &status);
        mpi_err(__FUNCTION__, ": After Waitany: ", pendingRequests);
        pendingRequests.pop_back();

        if (index != pendingRequests.size())
        {
            popped_special = false;
            return {status, popAt(index)};
        }
        else
        {
            popped_special = true;
            return {status, special_tag};
        }
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

    ~TaggedRequestQueue() { cancel_remaining(); }
};
