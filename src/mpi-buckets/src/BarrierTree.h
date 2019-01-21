#include "DataChannel.h"
#include "RequestQueue.h"

class BarrierTree
{
    DataChannel  channel;
    RequestQueue notifyQueue;
    int          my_rank;
    int          num_ranks;
    int          unfinished_children;

    bool finishedFlag = false;

    void notifyParents()
    {
        auto NotifyRank
            = [](DataChannel const& channel, int my_rank, RequestQueue& queue) {
                  queue.push(channel.notify(my_rank));
              };
        onParentRanks(NotifyRank, notifyQueue);
    }
    void markComplete() { markChildComplete(); }
    void notifyChildren()
    {
        auto NotifyRank
            = [](DataChannel const& channel, int my_rank, RequestQueue& queue) {
                  queue.push(channel.notify(my_rank));
              };
        finishedFlag = true;
        onChildRanks(NotifyRank, notifyQueue);
    }
    
    void markChildComplete()
    {
        unfinished_children -= 1;
        if (unfinished_children == 0)
        {
            if (isRoot())
            {
                notifyChildren();
            }
            else
            {
                notifyParents();
            }
        }
    }
    template <class F, class... ExtraArgs>
    auto onChildRanks(F&& func, ExtraArgs&&... args) const -> void
    {
        auto const child_0 = (my_rank + 1) * 2 - 1;
        auto const child_1 = (my_rank + 1) * 2;
        if (child_0 < num_ranks)
        {
            func(channel, (int)child_0, args...);
        }
        if (child_1 < num_ranks)
        {
            func(channel, (int)child_1, args...);
        }
    }
    template <class F, class... ExtraArgs>
    auto onParentRanks(F&& func, ExtraArgs&&... args) const -> void
    {
        if (my_rank != 0)
        {
            const int parent = (my_rank - 1) / 2;
            func(channel, parent, args...);
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
        int const parent = (my_rank - 1) / 2;
        return sender == parent;
    }
    auto isFromChild(int sender) const -> bool
    {
        return (sender - 1) / 2 == my_rank;
    }
    auto isRoot() const -> bool { return my_rank == 0; }
    auto countChildren() const -> int
    {
        int count = 0;
        onChildRanks([&](DataChannel const&, int) { ++count; });
        return count;
    }
    auto tag() const -> int { return channel.tag; }
    auto test() const -> bool { return finishedFlag; }

   public:
    BarrierTree(MPI_Comm comm, int tag)
      : channel({comm, tag})
      , notifyQueue()
      , my_rank(channel.commRank())
      , num_ranks(channel.commSize())
      // Here, we count ourselves as an unfinished child
      , unfinished_children(countChildren() + 1)
    {
    }    
    
    void processStatus(MPI_Status const& status)
    {
        if (isFromChild(status.MPI_SOURCE))
        {
            markChildComplete();
        }
        else if (isFromParent(status.MPI_SOURCE))
        {
            notifyChildren();
        }
    }
    void enqueueBarrier(RequestQueue& queue)
    {
        auto EnqueueAwait
            = [](DataChannel const& channel, int rank, RequestQueue& queue) {
                  queue.push(channel.Iawait(rank));
              };
        onParentRanks(EnqueueAwait, queue);
        onParentRanks(EnqueueAwait, queue);
    }
};
