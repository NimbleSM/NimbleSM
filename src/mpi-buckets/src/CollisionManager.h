#pragma once
#include "BarrierTree.h"
#include "BoundingBox.h"
#include "DataChannel.h"
#include "GridCell.h"
#include "GridHash.h"
#include "GridIndex.h"
#include "RequestQueue.h"
#include "WaitAnyResult.h"
#include "meta.h"

#include <unordered_map>
#include <unordered_set>
#include <utility>

#include "mpi_err.h"

// Takes a list of points from a kokkos view and constructs a map of
// grid indicies and the corresponding bounding box
template <class View, class Map = std::unordered_map<GridIndex, BoundingBox>>
auto gatherIntoBoundingBoxes(View&& values, double const cell_size) -> Map
{
    Map          boundingBoxLookup;
    size_t const count = values.extent(0) / 3;
    if (count > 0)
    {
        double const scale = 1.0 / cell_size;

        Point3d point{values(0), values(1), values(2)};

        auto current_index       = GridIndex(point, scale);
        auto current_cell_bounds = GridCellBounds(current_index, cell_size);
        auto current_box         = BoundingBox(point);

        for (size_t i = 1; i < count; i++)
        {
            point = {values(3 * i), values(3 * i + 1), values(3 * i + 2)};

            if (current_cell_bounds.contains(point, cell_size))
            {
                current_box.include(point);
            }
            else
            {
                // We include the current bounding box in the map,
                // then calculate a new bounding box, index, and cell
                current_box.includeInMap(boundingBoxLookup, current_index);
                current_box.centerAt(point);
                current_index = GridIndex(point, scale);
                current_cell_bounds.moveTo(current_index, cell_size);
            }
        }
        // Ensure that the most recent version of the current box is included in
        // the map
        current_box.includeInMap(boundingBoxLookup, current_index);
    }
    return boundingBoxLookup;
}


template <class PacketMap>
void SyncSendMessageSizes(PacketMap&&   packets,
                          DataChannel   channel,
                          RequestQueue& queue,
                          size_t*       sizeStorage)
{
    for (auto&& packet : packets)
    {
        auto  destination = packet.first;
        auto& message     = packet.second;

        *sizeStorage = message.size();
        queue.push(channel.Issend(sizeStorage, 1, destination));
        sizeStorage += 1;
    }
}

template <class PacketMap, class SourceIdentifiedCallback, class SendProcessedCallback>
auto IdentifySources(PacketMap&&              packets,
                     MPI_Comm                 comm,
                     int                      dataTag,
                     int                      barrierTag,
                     SourceIdentifiedCallback notifySourceIdentified,
                     SendProcessedCallback    notifySendProcessed) -> void
{
    DataChannel channel{comm, dataTag};
    DataChannel barrierChannel{comm, barrierTag};
    using namespace std;

    auto queue                 = RequestQueue();
    auto messageCounts         = vector<size_t>(packets.size());
    auto active_outgoing_count = size_t(messageCounts.size());
    auto barrier               = BarrierTree(barrierChannel);
    SyncSendMessageSizes(packets, channel, queue, messageCounts.data());
    barrier.enqueueBarrier(queue);

    auto incoming_size       = size_t();
    auto active_recv_request = channel.Irecv(&incoming_size, 1, MPI_ANY_SOURCE);
    queue.push(active_recv_request);

    if (active_outgoing_count == 0)
    {
        barrier.markComplete();
    }
    mpi_err("Entering barrier loop; queue has ", queue.size(), " items");

    while (not barrier.test())
    {
        mpi_err("queue.pop()");
        auto reply = queue.pop();
        if (reply.request == active_recv_request)
        {
            mpi_err("~queue.pop() active_recv_request: reply from ",
                    reply.status.MPI_TAG);
            notifySourceIdentified(reply.status.MPI_SOURCE, incoming_size);

            active_recv_request = channel.Irecv(&incoming_size, 1, MPI_ANY_SOURCE);
            queue.push(active_recv_request);
        }
        else if (reply.status.MPI_TAG == channel.tag)
        {
            mpi_err("~queue.pop() send processed");
            notifySendProcessed(reply.status.MPI_SOURCE);
            active_outgoing_count--;
            if (active_outgoing_count == 0)
            {
                barrier.markComplete();
            }
        }
        else if (reply.status.MPI_TAG == barrier.tag())
        {
            mpi_err("~queue.pop(): barrier.processStatus");
            barrier.processStatus(reply.status);
        }
        else
        {
            mpi_err("~queue.pop(): unknown");
        }
    }

    mpi_err("barrier complete");
}

template <class PacketMap,
          class DataT = ElemType<decltype(std::begin(iof<PacketMap&>())->second)>>
auto ExchangeData(PacketMap&& packets, MPI_Comm comm, int tag1, int tag2, int tag3)
    -> std::vector<std::pair<int, std::vector<DataT>>>
{
    using namespace std;

    auto incomingData = vector<pair<int, vector<DataT>>>();
    incomingData.reserve(10);

    auto         messageChannel = DataChannel({comm, tag3});
    RequestQueue queue;

    // clang-format off
    IdentifySources(
        packets,
        comm, 
        tag1, 
        tag2,
        [&](int source, size_t incoming_size) {
            incomingData.emplace_back();
            auto& back   = incomingData.back();
            back.first   = source;
            auto& buffer = back.second;
            buffer.resize(incoming_size);
            queue.push(messageChannel.Irecv(buffer.data(), incoming_size, source));
        },
        [&](int dest) {
            auto& message = packets[dest];
            queue.push(messageChannel.Isend(message, dest));
        }
    );
    mpi_err("Finished IdentifySources()");
    // clang-format on
    queue.wait_all();

    return incomingData;
}

auto anyIntersect(std::vector<BoundingBox> const& boxes1,
                  std::vector<BoundingBox> const& boxes2) -> bool
{
    for (auto b1 : boxes1)
    {
        for (auto const& b2 : boxes2)
        {
            if (b1.intersects(b2))
                return true;
        }
    }
    return false;
}

using BBPacket     = std::pair<int, std::vector<BoundingBox>>;
using BBPacketList = std::vector<BBPacket>;

/**
 * rankBoxInfo: list of all the bounding boxes my rank is responsible for,
 * as well as the ranks those bounding boxes came from
 *
 * sourceRanks: list of ranks from which my rank expects to recieve collision
 * info back from in response
 */
auto notifyRanksOfIntersection(BBPacketList const&     rankBoxInfo,
                               MPI_Comm                comm,
                               int                     tag1,
                               int                     tag2,
                               std::vector<int> const& sourceRanks)
    -> std::vector<int>
{
    auto dataChannel = DataChannel{comm, tag1};
    auto sizeChannel = DataChannel{comm, tag2};
    using namespace std;
    RequestQueue                 sizeRecvQueue;
    std::vector<size_t>          incomingSizes(sourceRanks.size());
    std::unordered_map<int, int> inverseSourceRanks = invert(sourceRanks);
    vector<vector<int>>          incomingRanks(sourceRanks.size());

    inverseSourceRanks.reserve(sourceRanks.size());

    // Get outgoingSizes of incoming messages
    for (int i = 0; i < incomingSizes.size(); i++)
    {
        sizeRecvQueue.push(sizeChannel.Irecv(&incomingSizes[i], 1, sourceRanks[i]));
    }

    RequestQueue sendQueue;
    auto         outgoingSizes     = vector<size_t>(rankBoxInfo.size());
    auto         intersectingRanks = vector<vector<int>>(rankBoxInfo.size());


    for (int index = 0; index < rankBoxInfo.size(); index++)
    {
        auto& destInfo    = rankBoxInfo[index];
        int   destRank    = destInfo.first;
        auto& listOfRanks = intersectingRanks[index];

        for (auto& otherInfo : rankBoxInfo)
        {
            if (&destInfo == &otherInfo)
                continue;
            if (anyIntersect(destInfo.second, otherInfo.second))
                listOfRanks.push_back(otherInfo.first);
        }

        sendQueue.push(sizeChannel.Isend(&outgoingSizes[index], 1, destRank));
        sendQueue.push(dataChannel.Isend(listOfRanks, destRank));
    }


    RequestQueue rankRecvQueue;
    size_t       total_incoming_ranks = 0;
    while (sizeRecvQueue.has())
    {
        WaitAnyResult result = sizeRecvQueue.pop();
        if (result.status.MPI_TAG == sizeChannel.tag)
        {
            int   source       = result.status.MPI_SOURCE;
            int   source_index = inverseSourceRanks.at(source);
            auto& vect         = incomingRanks[source_index];
            vect.resize(incomingSizes[source_index]);
            rankRecvQueue.push(dataChannel.Irecv(vect.data(), vect.size(), source));
            total_incoming_ranks += vect.size();
        }
    }
    rankRecvQueue.wait_all();


    std::unordered_set<int> ranks;
    ranks.reserve(total_incoming_ranks);
    for (auto& rankList : incomingRanks)
    {
        ranks.insert(rankList.begin(), rankList.end());
    }
    return std::vector<int>(ranks.begin(), ranks.end());
}

template <class View>
auto getExchangeMembers(View&&       kokkos_view,
                        double const cell_size,
                        MPI_Comm     comm,
                        int          tag1,
                        int          tag2,
                        int          tag3) -> std::vector<int>
{
    using namespace std;
    // Steps:
    // 1) Get distribution of points as bounding boxes within grid cells
    // 2) Group bounding boxes based on which rank they should be sent to
    // 3 outgoing) Distribute my bounding boxes to the corresponding ranks
    // 3 incoming) Recieve incoming bounding boxes from other ranks
    // 4 outgoing)
    //   - For each bounding,

    HashFunction hash;
    DataChannel  channel = DataChannel{comm, tag1};
    int const    n_ranks = channel.commSize();

    // clang-format off

    // Creates a map from ranks to the bounding boxes that
    // have to be sent to those ranks
    auto packets = GatherTransformBy(
        gatherIntoBoundingBoxes(
            kokkos_view, 
            cell_size
        ),
        [=](pair<const GridIndex, BoundingBox> const& message) {
            auto index = message.first;
            return hash(index.x, index.y, index.z) % n_ranks;
        },
        [](pair<const GridIndex, BoundingBox> const& message) {
            return message.second; 
        }
    ); 
    auto rankBoxInfo =
        ExchangeData(
            packets,
            channel.comm,
            tag1, 
            tag2,
            tag3
        );

    // clang-format on

    // source_ranks will contain the ranks that my rank sent bounding box info
    // to; my rank expects to recieve a response containing the list of ranks my
    // rank intersects with
    std::vector<int> source_ranks;
    for (auto& packet : packets)
    {
        source_ranks.push_back(packet.first);
    }
    return notifyRanksOfIntersection(rankBoxInfo, comm, tag1, tag2, source_ranks);
}
