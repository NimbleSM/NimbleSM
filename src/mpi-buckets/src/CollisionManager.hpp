#pragma once
#include "BarrierTree.hpp"
#include "GatherIntoBoundingBoxes.hpp"

#include "GridIndex.hpp"
#include "RequestQueue.hpp"
#include "TaggedRequestQueue.hpp"
#include "WaitAnyResult.hpp"
#include "meta.hpp"

#include <unordered_map>
#include <unordered_set>
#include <utility>

#include "mpi_err.hpp"


enum class RequestType { Incoming, Outgoing, Barrier };



/** Given a bunch of outgoing packets, SyncSendMessageSizes sends out the
message sizes of the packets to the corresponding destination. It uses
channel.Issend because the request is meant to delay until the reciever
accepts the incoming size.*/
template <class PacketMap>
void SyncSendMessageSizes(PacketMap&&                      packets,
                          DataChannel                      channel,
                          TaggedRequestQueue<RequestType>& queue,
                          size_t*                          sizeStorage)
{
    for (auto&& packet : packets)
    {
        auto  destination = packet.first;
        auto& message     = packet.second;

        *sizeStorage = message.size();
        queue.push(channel.Issend(sizeStorage, 1, destination),
                   RequestType::Outgoing);
        sizeStorage += 1;
    }
}

/**
IdentifySources
*/
template <class PacketMap, class SourceIdentifiedCallback, class SendProcessedCallback>
auto IdentifySources(PacketMap&&              packets,
                     MPI_Comm                 comm,
                     int                      dataTag,
                     int                      barrierTag,
                     SourceIdentifiedCallback notifySourceIdentified,
                     SendProcessedCallback    notifySendProcessed) -> void
{
    using namespace std;

    DataChannel                     channel{comm, dataTag};
    BarrierTree                     barrier{DataChannel{comm, barrierTag}};
    TaggedRequestQueue<RequestType> queue{};
    vector<size_t>                  messageCounts(packets.size());
    size_t active_outgoing_count = messageCounts.size();

    barrier.enqueueBarrier(queue, RequestType::Barrier);
    if (active_outgoing_count == 0)
    {
        barrier.markComplete();
    }

    SyncSendMessageSizes(packets, channel, queue, messageCounts.data());

    size_t incoming_size     = 0;
    auto active_recv_request = channel.Irecv(&incoming_size, 1, MPI_ANY_SOURCE);
    int  active_recieved     = 0;

    while (not barrier.test())
    {
        mpi_err(__FUNCTION__, ": items in queue: ", queue.pendingRequests);
        bool active_request_fulfilled = false;
        mpi_err("active request: ", active_recv_request);
        auto reply = queue.popWithSpecial(
            active_recv_request, RequestType::Incoming, active_request_fulfilled);
        mpi_err("active request: ", active_recv_request);
        if (active_request_fulfilled)
        {
            mpi_err("~queue.pop() active_recv_request: reply from ",
                    reply.first.MPI_TAG);
            notifySourceIdentified(reply.first.MPI_SOURCE, incoming_size);

            active_recv_request = channel.Irecv(&incoming_size, 1, MPI_ANY_SOURCE);
        }
        else if (reply.second == RequestType::Outgoing)
        {
            mpi_err("~queue.pop() send processed from ", reply.first.MPI_SOURCE);
            notifySendProcessed(reply.first.MPI_SOURCE);
            active_outgoing_count--;
            if (active_outgoing_count == 0)
            {
                barrier.markComplete();
            }
        }
        else if (reply.second == RequestType::Barrier)
        {
            mpi_err("~queue.pop(): barrier.processStatus");
            barrier.processStatus(reply.first);
        }
        else
        {
            mpi_err("~queue.pop(): unknown");
        }
    }
    mpi_err("barrier complete");
    {
        MPI_Status reply{};
        int        flag = 0;
        MPI_Test(&active_recv_request, &flag, &reply);
        if (flag)
        {
            mpi_err("~queue.pop() active_recv_request: reply from ", reply.MPI_TAG);
            notifySourceIdentified(reply.MPI_SOURCE, incoming_size);
        }
    }
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
            mpi_err("ExchangeData recv from ", source); 
            queue.push(messageChannel.Irecv(buffer.data(), incoming_size, source));
        },
        [&](int dest) {
            mpi_err("ExchangeData send to ", dest); 
            auto& message = packets[dest];
            queue.push(messageChannel.Isend(message, dest));
        }
    );
    mpi_err("Finished IdentifySources()");
    // clang-format on
    queue.wait_all();
    mpi_err("finished waiting on queue");

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
    TaggedRequestQueue<RequestType> sizeRecvQueue;

    std::vector<size_t>          incomingSizes(sourceRanks.size());
    std::unordered_map<int, int> inverseSourceRanks = invert(sourceRanks);
    vector<vector<int>>          incomingRanks(sourceRanks.size());

    inverseSourceRanks.reserve(sourceRanks.size());

    // Get outgoingSizes of incoming messages
    mpi_err("notifyRanksOfIntersection(): sourceRanks: ", sourceRanks);
    mpi_err("notifyRanksOfintersection(): rankBoxInfo.size(): ",
            rankBoxInfo.size());
    for (int i = 0; i < incomingSizes.size(); i++)
    {
        sizeRecvQueue.push(sizeChannel.Irecv(&incomingSizes[i], 1, sourceRanks[i]),
                           RequestType::Incoming);
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

        outgoingSizes[index] = listOfRanks.size();
        mpi_err("Sending stuff to ", destRank);
        sendQueue.push(sizeChannel.Isend(&outgoingSizes[index], 1, destRank));
        sendQueue.push(dataChannel.Isend(listOfRanks, destRank));
    }


    RequestQueue rankRecvQueue;
    size_t       total_incoming_ranks = 0;
    mpi_err("notifyRanksOfIntersection(): entering sizeRecvQueue loop");
    while (sizeRecvQueue.has())
    {
        mpi_err("sizeRecvQueue.pop()");
        auto       tag_result_pair = sizeRecvQueue.pop();
        MPI_Status result          = tag_result_pair.first;
        auto       category        = tag_result_pair.second;
        mpi_err("~sizeRecvQueue.pop()");
        if (result.MPI_TAG == sizeChannel.tag)
        {
            int   source       = result.MPI_SOURCE;
            int   source_index = inverseSourceRanks.at(source);
            auto& vect         = incomingRanks[source_index];
            vect.resize(incomingSizes[source_index]);
            rankRecvQueue.push(dataChannel.Irecv(vect.data(), vect.size(), source));
            total_incoming_ranks += vect.size();
        }
    }
    mpi_err("notifyRanksOfIntersection(): rankRecvQueue.wait_all()");
    rankRecvQueue.wait_all();
    mpi_err("~rankRecvQueue.wait_all()");

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
    mpi_err("ExchangeData()"); 
    auto rankBoxInfo =
        ExchangeData(
            packets,
            channel.comm,
            tag1, 
            tag2,
            tag3
        );
    
    mpi_err("~ExchangeData()");

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
