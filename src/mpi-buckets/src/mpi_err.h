#pragma once
#include <mpi.h>
#include <cstdio>
#include <string>


template <size_t N>
void str(std::string& s, char const (&msg)[N])
{
    s.append(msg, N);
}
void str(std::string& s, const char* msg) { s.append(msg); }
void str(std::string& s, short i) { s.append(std::to_string(i)); }
void str(std::string& s, int i) { s.append(std::to_string(i)); }
void str(std::string& s, long int i) { s.append(std::to_string(i)); }
void str(std::string& s, unsigned i) { s.append(std::to_string(i)); }
void str(std::string& s, unsigned long i) { s.append(std::to_string(i)); }
void str(std::string& s, bool b)
{
    if (b)
        s.append("true", 4);
    else
        s.append("false", 5);
}

int getMPIRank()
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    return rank;
}
template <class... T>
void mpi_err(T&&... things)
{
    // static std::string rank = "rank " + std::to_string(getMPIRank()) + ": ";
    // std::string        s    = rank;
    // char               _[]{(str(s, things), '\0')...};
    // s.push_back('\n');
    // fwrite(s.data(), 1, s.size(), stdout);
    // fflush(stderr);
}
