#pragma once
#include <GatherIntoBoundingBoxes.hpp>
#include <mpi.h>
#include <mpi/DataChannel.hpp>
#include <vector>

template<class T> struct simple_span
{
    T* first;
    T* last;
    simple_span(T* f, T* l)
      : first(f)
      , last(l)
    {}
    simple_span(T* f, size_t count)
      : first(f)
      , last(f + count)
    {}

    size_t size() const { return last - first; }
    T* begin() const { return first; }
    T* end() const { return last; }
};
template<class T> simple_span<T> make_span(T* ptr, size_t size)
{
    return simple_span<T>(ptr, size);
}
struct grid_bb_pair
{
    GridIndex index;
    BoundingBox box;
};
template<class View>
auto get_exchange_members_root(View&& kokkos_view,
                               double const cell_size,
                               MPI_Comm comm,
                               int tag1,
                               int tag2,
                               int tag3) -> std::vector<int>
{
    std::unordered_map<GridIndex, BoundingBox> map =
      gather_into_bounding_boxes(kokkos_view, cell_size);


    std::vector<grid_bb_pair> map_buffer(map.size());

    {
        auto* map_scan = map_buffer.data();
        for (auto& elem : map)
        {
            map_scan->index = elem.first;
            map_scan->box   = elem.second;
            map_scan++;
        }
    }

    DataChannel channel {comm, tag1};

    int my_rank             = channel.commRank();
    int num_ranks           = channel.commSize();
    int outgoing_size_bytes = map_buffer.size() * sizeof(grid_bb_pair);
    if (my_rank == 0)
    {
        std::vector<int> incoming_sizes_bytes(num_ranks);
        std::vector<int> incoming_counts(num_ranks);

        MPI_Gather(&outgoing_size_bytes,
                   1,
                   MPI_INT,
                   incoming_sizes_bytes.data(),
                   1,
                   MPI_INT,
                   0,
                   comm);

        int total_incoming = 0;
        std::vector<int> incoming_displacements_bytes(num_ranks);
        std::vector<int> incoming_displacements(num_ranks);

        {
            int* disp_bytes = incoming_displacements_bytes.data();
            int* disp_count = incoming_displacements.data();
            int* count_ptr  = incoming_counts.data();
            for (int size_bytes : incoming_sizes_bytes)
            {
                *disp_bytes = total_incoming;
                *disp_count = total_incoming / sizeof(grid_bb_pair);
                *count_ptr  = size_bytes / sizeof(grid_bb_pair);
                count_ptr++;
                disp_count++;
                disp_bytes++;

                total_incoming += size_bytes;
            }
        }
        std::vector<grid_bb_pair> incoming_buffer(total_incoming /
                                                  sizeof(grid_bb_pair));

        MPI_Gatherv(map_buffer.data(),
                    outgoing_size_bytes,
                    MPI_BYTE,
                    incoming_buffer.data(),
                    incoming_sizes_bytes.data(),
                    incoming_displacements_bytes.data(),
                    MPI_BYTE,
                    0,
                    comm);

        std::vector<std::vector<int>> intersecting(num_ranks);

        for (int i = 0; i < num_ranks; i++)
        {
            auto s1 =
              make_span(incoming_buffer.data() + incoming_displacements[i],
                        incoming_counts[i]);

            for (int j = 0; j < num_ranks; j++)
            {
                auto s2 =
                  make_span(incoming_buffer.data() + incoming_displacements[j],
                            incoming_counts[i]);

                for (auto& elem1 : s1)
                {
                    for (auto& elem2 : s2)
                    {
                        if (elem1.index == elem2.index &&
                            elem1.box.intersects(elem2.box))
                        {
                            intersecting[i].push_back(j);
                            goto continue_outer;
                        }
                    }
                }

            continue_outer:
                continue;
            }
        }

        // Scatter the ranks back out
    }
    else
    {
        MPI_Gather(
          &outgoing_size_bytes, 1, MPI_INT, nullptr, 0, MPI_INT, 0, comm);

        MPI_Gatherv(map_buffer.data(),
                    outgoing_size_bytes,
                    MPI_BYTE,
                    nullptr,
                    nullptr,
                    nullptr,
                    MPI_BYTE,
                    0,
                    comm);
    }
    return {};
}

template<class View>
auto get_exchange_members_basic(View&& kokkos_view,
                                double const cell_size,
                                MPI_Comm comm,
                                int tag1,
                                int tag2,
                                int tag3) -> std::vector<int>
{
    DataChannel channel = DataChannel {comm, tag1};
    int const n_ranks   = channel.commSize();

    std::vector<int> ranks(n_ranks);
    for (int i = 0; i < n_ranks; i++)
    {
        ranks[i] = i;
    }

    return ranks;
}

template<class View>
auto get_echange_members(View&& kokkos_view,
                         double const cell_size,
                         MPI_Comm comm,
                         int t1,
                         int t2,
                         int t3) -> std::vector<int>
{
    return get_exchange_members_basic(kokkos_view, cell_size, comm, t1, t2, t3);
}