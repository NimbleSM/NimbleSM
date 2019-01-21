#include "DataChannel.h"
#include "RequestQueue.h"

class BarrierTree
{
    int          unfinished_children;
    int          rank;
    DataChannel  channel;
    RequestQueue notifyQueue;

   public:
    void signalComplete()
    {
        constexpr auto NotifyRank
            = [](DataChannel const& channel, int rank, RequestQueue& queue) {
                  queue.push(channel.notify(rank));
              };
        if (unfinished_children == 0)
        {
            if (channel.isRoot())
            {
                channel.onChildRanks(NotifyRank, notifyQueue);
            }
            else
            {
                channel.onParentRanks(NotifyRank, notifyQueue);
            }
        }
    }
}; 
