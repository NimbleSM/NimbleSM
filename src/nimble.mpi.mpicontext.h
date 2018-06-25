/*
//@HEADER
// ************************************************************************
//
//                                NimbleSM
//                             Copyright 2018
//   National Technology & Engineering Solutions of Sandia, LLC (NTESS)
//
// Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government
// retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY EXPRESS OR IMPLIED
// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
// MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN
// NO EVENT SHALL NTESS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
// TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions?  Contact David Littlewood (djlittl@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#pragma once
#include <mpi.h>

#include <algorithm>
#include <array>
#include <cstring>
#include <exception>
#include <iostream>
#include <memory>
#include <numeric>
#include <string>
#include <tuple>
#include <type_traits>
#include <vector>

#include "nimble.mpi.serialization.h"

namespace nimble
{
class mpicontext
{
  static std::vector<int> serialization_buffer;
  static int get_rank(MPI_Comm comm)
  {
    int rank;
    MPI_Comm_rank(comm, &rank);
    return rank;
  }
  static int get_size(MPI_Comm comm)
  {
    int size;
    MPI_Comm_size(comm, &size);
    return size;
  }
  int rank = 0;
  int size = 0;
  MPI_Comm comm;

 public:
  const MPI_Comm& get_comm() const { return comm; }
  int get_size() const { return size; }
  int get_rank() const { return rank; }
  mpicontext() = default;
  mpicontext(MPI_Comm comm) : rank{get_rank(comm)}, size{get_size(comm)}, comm{comm} {}
  mpicontext(const mpicontext& context) = default;
  mpicontext(mpicontext&& context)      = default;
  mpicontext& operator=(const mpicontext& other) = default;
  mpicontext& operator=(mpicontext&& other) = default;
  int get_root() const { return 0; }
  bool is_root() const { return rank == get_root(); }
  void send(const std::vector<int>& data, int rank, int tag) const
  {
    MPI_Send(data.data(), data.size(), MPI_INT, rank, tag, this->comm);
  }
  void send(std::pair<int*, int> data, int rank, int tag) const
  {
    MPI_Send(data.first, data.second, MPI_INT, rank, tag, this->comm);
  }
  int recv_count(int rank, int tag) const
  {
    int count = 0;
    MPI_Status recv_status;
    MPI_Probe(rank, tag, this->comm, &recv_status);
    MPI_Get_count(&recv_status, MPI_INT, &count);
    return count;
  }
  int recv(std::vector<int>& data, int rank, int tag) const
  {
    int count = recv_count(rank, tag);
    data.resize(count);
    MPI_Recv(data.data(), count, MPI_INT, rank, tag, this->comm, MPI_STATUS_IGNORE);
    return count;
  }
  int recv(int* data, int rank, int tag, int max_buffer_size) const
  {
    MPI_Status recv_status;
    MPI_Recv(data, max_buffer_size, MPI_INT, rank, tag, this->comm, &recv_status);
    int count;
    MPI_Get_count(&recv_status, MPI_INT, &count);
    return count;
  }
  // This will only resize the vector if it needs more space, but not if it's smaller
  int recv_avoid_resize(std::vector<int>& data, int rank, int tag) const
  {
    int count = recv_count(rank, tag);
    if (count > data.size())
      data.resize(count);
    MPI_Recv(data.data(), count, MPI_INT, rank, tag, this->comm, MPI_STATUS_IGNORE);
    return count;
  }
  template<class T, class podT>
  void sendpacked(const T& source, int rank, int tag, std::vector<podT>& buffer) const
  {
    serialization::pack_avoid_resize(source, buffer);
    send(buffer, rank, tag);
  }
  template<class T, class podT>
  void recvpacked(T& dest, int rank, int tag, std::vector<podT>& buffer) const
  {
    recv_avoid_resize(buffer, rank, tag);
    serialization::unpack(dest, buffer);
  }

  // Attempts to serialize the input and send it.
  template<class T>
  void sendpacked(const T& source, int rank, int tag) const
  {
    int sendcount = serialization::pack_avoid_resize(source, serialization_buffer);
    send({serialization_buffer.data(), sendcount}, rank, tag);
  }
  std::vector<std::vector<int>> gather(const std::vector<int>& source, int root) const
  {
    if (this->get_rank() == root)
    {
      // Step 1: figure out how much data is going to be incoming from each rank
      std::vector<int> recv_counts(this->get_size());
      int sendcount = source.size();
      // TO-DO: review documentation at
      // https://www.mpich.org/static/docs/v3.1/www3/MPI_Gather.html
      // https://www.mpich.org/static/docs/v3.1/www3/MPI_Gatherv.html
      // to ensure that the function call is being passed the correct aruments.
      MPI_Gather(&sendcount, 1, MPI_INT, recv_counts.data(), 1, MPI_INT, root, this->comm);

      // Step 2: Get data from all ranks as contiguous block of data
      // destination.size() is the sum of the lengths of all incoming arrays.
      std::vector<int> destination(std::accumulate(recv_counts.begin(), recv_counts.end(), 0));
      std::vector<int> displacements(this->get_size(), 0);
      std::partial_sum(recv_counts.begin(), recv_counts.end() - 1, displacements.data() + 1);
      MPI_Gatherv(source.data(),
                  sendcount,
                  MPI_INT,
                  destination.data(),
                  recv_counts.data(),
                  displacements.data(),
                  MPI_INT,
                  root,
                  this->get_comm());

      // Step 3: Copy data from a contiguous block into a vector of vectors.
      std::vector<std::vector<int>> result(this->get_size());
      int *start_iter = destination.data(), *end_iter = destination.data() + recv_counts[0];
      for (size_t i = 0; i < result.size(); ++i)
      {
        result[i]  = std::vector<int>(start_iter, end_iter);
        start_iter = end_iter;
        end_iter += recv_counts[i];
      }
      // Return result
      return result;
    }
    else
    {
      // Step 1 (see above)
      int sendcount = source.size();
      MPI_Gather(&sendcount, 1, MPI_INT, 0, 0, MPI_INT, root, this->get_comm());

      // Step 2 (see above)
      MPI_Gatherv(source.data(), sendcount, MPI_INT, 0, 0, 0, MPI_INT, root, this->get_comm());

      // Return empty vector, as only the root rank will recieve incoming data.
      return {};
    }
  }
  std::vector<int> scatter(int root, const std::vector<std::vector<int>>& source) const
  {
    if (this->get_rank() == root)
    {
      // Step 1: notify each rank of how much data they should expect to recieve
      int recvcount;
      std::vector<int> send_counts(this->get_size());
      for (size_t i = 0; i < send_counts.size(); ++i)
      {
        send_counts[i] = source[i].size();
      }
      MPI_Scatter(send_counts.data(), 1, MPI_INT, &recvcount, 1, MPI_INT, root, this->get_comm());

      // Step 2: send out data to each of the ranks.
      std::vector<int> destination(recvcount, 0);
      std::vector<int> displacements(this->get_size(), 0);
      std::vector<int> source_made_contiguous(
          std::accumulate(send_counts.begin(), send_counts.end(), 0), 0);

      std::partial_sum(send_counts.begin(), send_counts.end() - 1, displacements.data() + 1);
      MPI_Scatterv(source.data(),
                   send_counts.data(),
                   displacements.data(),
                   MPI_INT,
                   destination.data(),
                   recvcount,
                   MPI_INT,
                   root,
                   this->get_comm());
    }
    else
    {
      int recvcount;
      MPI_Scatter(0, 0, MPI_INT, &recvcount, 1, MPI_INT, root, this->get_comm());
      std::vector<int> recvbuff(recvcount);
      MPI_Scatterv(0, 0, 0, MPI_INT, recvbuff.data(), recvcount, MPI_INT, root, this->get_comm());
      return recvbuff;
    }
    return std::vector<int>();
  }

  template<size_t count>
  std::vector<std::array<int, count>> allgather_to_vector(
      const std::array<int, count>& values) const
  {
    std::vector<std::array<int, count>> buffer(get_size());
    MPI_Allgather(values.data(),      /* data being sent */
                  (int)values.size(), /* amount of data being sent */
                  MPI_INT,            /* Send type */
                  buffer.data(),      /* buffer to place revieced data */
                  (int)values.size(), /* spacing with which to insert data from buffer */
                  MPI_INT,            /* Recieve type */
                  get_comm()          /* MPI Comm to use */
    );
    return buffer;
  }
  std::vector<int> allgather_to_vector(int value) const
  {
    std::vector<int> buffer(get_size());
    MPI_Allgather(&value,        /* data being sent */
                  1,             /* amount of data being sent */
                  MPI_INT,       /* Send type */
                  buffer.data(), /* buffer to place revieced data */
                  1,             /* spacing with which to insert data from buffer */
                  MPI_INT,       /* Recieve type */
                  get_comm()     /* MPI Comm to use */
    );
    return buffer;
  }
  std::vector<size_t> allgather_to_vector(size_t value) const
  {
    std::vector<size_t> buffer(get_size());
    if (sizeof(unsigned long) == sizeof(size_t))
    {
      MPI_Allgather(&value,            /* data being sent */
                    1,                 /* amount of data being sent */
                    MPI_UNSIGNED_LONG, /* Send type */
                    buffer.data(),     /* buffer to place revieced data */
                    1,                 /* spacing with which to insert data from buffer */
                    MPI_UNSIGNED_LONG, /* Recieve type */
                    get_comm()         /* MPI Comm to use */
      );
    }
    else
    {
      MPI_Allgather(&value,            /* data being sent */
                    sizeof(size_t),    /* amount of data being sent */
                    MPI_UNSIGNED_CHAR, /* Send type */
                    buffer.data(),     /* buffer to place revieced data */
                    sizeof(size_t),    /* spacing with which to insert data from buffer */
                    MPI_UNSIGNED_CHAR, /* Recieve type */
                    get_comm()         /* MPI Comm to use */
      );
    }
    return buffer;
  }
  void gather_send(int value) const
  {
    MPI_Gather(&value,     /* data being sent */
               1,          /* amount of data being sent */
               MPI_INT,    /* Send type */
               nullptr,    /* buffer to place revieced data */
               1,          /* spacing with which to insert data from buffer */
               MPI_INT,    /* Recieve type */
               get_root(), /* MPI Comm to use */
               get_comm()  /* rank that acts as root */
    );
  }
  std::vector<int> gather_recieve(int value) const
  {
    std::vector<int> buffer(get_size());
    MPI_Gather(&value,        /* data being sent */
               1,             /* amount of data being sent */
               MPI_INT,       /* Send type */
               buffer.data(), /* buffer to place revieced data */
               1,             /* spacing with which to insert data from buffer */
               MPI_INT,       /* Recieve type */
               get_root(),    /* MPI Comm to use */
               get_comm()     /* rank that acts as root */
    );
    return buffer;
  }
  template<size_t n>
  void gather_send(const std::array<int, n>& arr) const
  {
    MPI_Gather(arr.data(), /* data being sent */
               n,          /* amount of data being sent */
               MPI_INT,    /* Send type */
               nullptr,    /* buffer to place revieced data */
               n,          /* spacing with which to insert data from buffer */
               MPI_INT,    /* Recieve type */
               get_root(), /* MPI Comm to use */
               get_comm()  /* rank that acts as root */
    );
  }
  template<size_t n>
  std::vector<std::array<int, n>> gather_recieve(const std::array<int, n>& arr) const
  {
    std::vector<std::array<int, n>> buffer(get_size());
    MPI_Gather(arr.data(),    /* data being sent */
               n,             /* amount of data being sent */
               MPI_INT,       /* Send type */
               buffer.data(), /* buffer to place revieced data */
               n,             /* spacing with which to insert data from buffer */
               MPI_INT,       /* Recieve type */
               get_root(),    /* MPI Comm to use */
               get_comm()     /* rank that acts as root */
    );
    return buffer;
  }
  void gather_send(int* data, int size) const
  {
    MPI_Gather(data,       /* data being sent */
               size,       /* amount of data being sent */
               MPI_INT,    /* Send type */
               nullptr,    /* buffer to place revieced data */
               size,       /* spacing with which to insert data from buffer */
               MPI_INT,    /* Recieve type */
               get_root(), /* MPI Comm to use */
               get_comm()  /* rank that acts as root */
    );
  }
  std::vector<int> gather_recieve(int* data, int size) const
  {
    std::vector<int> buffer(get_size() * size);
    MPI_Gather(data,          /* data being sent */
               size,          /* amount of data being sent */
               MPI_INT,       /* Send type */
               buffer.data(), /* buffer to place revieced data */
               size,          /* spacing with which to insert data from buffer */
               MPI_INT,       /* Recieve type */
               get_root(),    /* MPI Comm to use */
               get_comm());   /* rank that acts as root */
    return buffer;
  }

  void bcast(int& value) const
  {
    MPI_Bcast(&value, 1, MPI_INT, this->get_root(), this->get_comm());
  }
  template<size_t count>
  void bcast(std::array<int, count>& arr) const
  {
    MPI_Bcast(arr.data(), count, MPI_INT, this->get_root(), this->get_comm());
  }

  std::string catenate_gatherv(const std::string& str) const
  {
    if (is_root())
    {
      std::vector<int> counts{gather_recieve(str.size())};
      std::vector<int> displacements(counts.size());
      std::partial_sum(counts.begin(), counts.end() - 1, displacements.begin() + 1);
      std::string dest((size_t)(displacements.back() + counts.back()), ' ');
      MPI_Gatherv(str.data(),           /* data being sent */
                  str.size(),           /* amount of data being sent */
                  MPI_CHAR,             /* type of data being sent */
                  &dest[0],             /* destination of data being sent */
                  counts.data(),        /* expected amount of data to be recieved from each rank */
                  displacements.data(), /* displacements of revieced data */
                  MPI_CHAR,             /* type of data being recieved */
                  this->get_root(),     /* root */
                  this->get_comm()      /* MPI Comm to use */
      );
      return dest;
    }
    else
    {
      gather_send(str.size());
      MPI_Gatherv(str.data(), /* data being sent */
                  str.size(), /* amount of data being sent */
                  MPI_CHAR,   /* type of data being sent */
                  nullptr,    /* destination of data being sent */
                  nullptr,    /* expected amount of data to be recieved from each rank */
                  nullptr,    /* displacements of revieced data */
                  MPI_CHAR,   /* type of data being recieved */
                  get_root(), /* root */
                  get_comm()  /* MPI Comm to use */
      );
      return {};
    }
  }
  std::string catenate_gatherv_format(const std::string& str,
                                      const std::string& _start,
                                      const std::string& separator,
                                      const std::string& _end) const
  {
    if (is_root())
    {
      std::vector<int> counts{gather_recieve(str.size())};
      std::vector<int> displacements;
      displacements.reserve(counts.size());
      int displacement_accumulator = _start.size();
      int separator_size           = separator.size();
      for (int count : counts)
      {
        displacements.emplace_back(displacement_accumulator);
        displacement_accumulator += separator_size + count;
      }
      std::string dest((size_t)(displacement_accumulator + _end.size() - separator_size), ' ');
      MPI_Gatherv(str.data(),           /* data being sent */
                  str.size(),           /* amount of data being sent */
                  MPI_CHAR,             /* type of data being sent */
                  &dest[0],             /* destination of data being sent */
                  counts.data(),        /* expected amount of data to be recieved from each rank */
                  displacements.data(), /* displacements of revieced data */
                  MPI_CHAR,             /* type of data being recieved */
                  this->get_root(),     /* root */
                  this->get_comm()      /* MPI Comm to use */
      );
      std::copy(_start.begin(), _start.end(), &dest[0]);
      for (int i = 0, max = counts.size() - 1; i < max; ++i)
      {
        int offset = displacements[i] + counts[i];
        std::copy(separator.begin(), separator.end(), &dest[offset]);
      }
      std::copy(_end.begin(), _end.end(), &dest[displacements.back() + counts.back()]);
      return dest;
    }
    else
    {
      gather_send(str.size());
      MPI_Gatherv(str.data(), /* data being sent */
                  str.size(), /* amount of data being sent */
                  MPI_CHAR,   /* type of data being sent */
                  nullptr,    /* destination of data being sent */
                  nullptr,    /* expected amount of data to be recieved from each rank */
                  nullptr,    /* displacements of revieced data */
                  MPI_CHAR,   /* type of data being recieved */
                  get_root(), /* root */
                  get_comm()  /* MPI Comm to use */
      );
      return {};
    }
  }
  const mpicontext& print(const std::string& s, std::ostream& os = std::cout) const
  {
    std::string str = catenate_gatherv(s);
    if (is_root())
    {
      os << str;
    }
    return *this;
  }
  const mpicontext& print_formatted(const std::string& s,
                                    const std::string& _start,
                                    const std::string& separator,
                                    const std::string& _end,
                                    std::ostream& os = std::cout) const
  {
    std::string str = catenate_gatherv_format(s, _start, separator, _end);
    if (is_root())
    {
      os << str;
    }
    return *this;
  }
  const mpicontext& println(const std::string& s, std::ostream& os = std::cout) const
  {
    std::string str = catenate_gatherv(s + "\n");
    if (is_root())
    {
      os << str;
    }
    return *this;
  }
  template<class T>
  const mpicontext& print_if_root(const T& s, std::ostream& os = std::cout) const
  {
    if (is_root())
    {
      os << s;
    }
    return *this;
  }
  template<class T>
  const mpicontext& println_if_root(const T& s, std::ostream& os = std::cout) const
  {
    if (is_root())
    {
      os << s << std::endl;
    }
    return *this;
  }
  void gatherv_send(const std::vector<int>& ints) const
  {
    // See: https://www.open-mpi.org/doc/v2.1/man3/MPI_Gatherv.3.php
    MPI_Gatherv(ints.data(), /* data being sent */
                ints.size(), /* amount of data being sent */
                MPI_INT,     /* type of data being sent */
                nullptr,     /* destination of data being sent */
                nullptr,     /* expected amount of data to be recieved from each rank */
                nullptr,     /* displacements of revieced data */
                MPI_INT,     /* type of data being recieved */
                get_root(),  /* root */
                get_comm()); /* MPI Comm to use */
  }
  void gatherv_recieve(const std::vector<int>& ints,
                       std::vector<int>& dest,
                       const std::vector<int>& counts,
                       const std::vector<int>& displacements) const
  {
    // See: https://www.open-mpi.org/doc/v2.1/man3/MPI_Gatherv.3.php
    MPI_Gatherv(ints.data(),          /* data being sent */
                ints.size(),          /* amount of data being sent */
                MPI_INT,              /* type of data being sent */
                dest.data(),          /* destination of data being sent */
                counts.data(),        /* expected amount of data to be recieved from each rank */
                displacements.data(), /* displacements of revieced data */
                MPI_INT,              /* type of data being recieved */
                this->get_root(),     /* root */
                this->get_comm());    /* MPI Comm to use */
  }

  void scatterv_send(const std::vector<int>& ints,
                     const std::vector<int>& counts,
                     const std::vector<int>& displacements,
                     std::vector<int>& dest) const
  {
    // See: https://www.open-mpi.org/doc/v2.1/man3/MPI_Scatterv.3.php
    MPI_Scatterv(ints.data(),          /* data being sent */
                 counts.data(),        /* amount of data being sent to each rank */
                 displacements.data(), /* displacements of data being sent */
                 MPI_INT,              /* type of data being sent */
                 dest.data(),          /* place to put sent data */
                 dest.size(),          /* amount of data to recieve */
                 MPI_INT,              /* type of data being recieved */
                 this->get_root(),     /* root sending the data */
                 this->get_comm());    /*communicator */
  }
  void scatterv_recieve(std::vector<int>& dest) const
  {
    // See: https://www.open-mpi.org/doc/v2.1/man3/MPI_Scatterv.3.php
    MPI_Scatterv(nullptr,           /* data being sent */
                 nullptr,           /* amount of data being sent to each rank */
                 nullptr,           /* displacements of data being sent */
                 MPI_INT,           /* type of data being sent */
                 dest.data(),       /* place to put sent data */
                 dest.size(),       /* amount of data to recieve */
                 MPI_INT,           /* type of data being recieved */
                 this->get_root(),  /* root sending the data */
                 this->get_comm()); /*communicator */
  }
  MPI_Comm split_by_color(int color) const { return split_by_color(color, this->get_rank()); }
  MPI_Comm split_by_color(int color, int key) const
  {
    MPI_Comm new_comm;
    MPI_Comm_split(get_comm(), color, key, &new_comm);
    return new_comm;
  }
  MPI_Comm split_by_color() const { return split_by_color(MPI_UNDEFINED); }
};
}   // namespace nimble