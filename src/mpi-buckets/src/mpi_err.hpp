#pragma once
#include <mpi.h>
#include <cstdio>
#include <string>
#include <vector>
#include <type_traits>


void str(std::string& s, void const* ptr)
{
    char buff[30];
    int  n = snprintf(buff, 30, "%p", ptr);
    s.append(buff, n);
}

template <size_t N>
void str(std::string& s, char const (&msg)[N])
{
    s.append(msg, N);
}
void str(std::string& s, const char* msg) { s.append(msg); }
void str(std::string& s, char* msg) { s.append(msg); }
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
template <class T>
void str(std::string& s, std::vector<T> const& vals)
{
    if (vals.size() == 0)
        s.append("[]");
    else
    {
        s.push_back('[');
        for (auto& elem : vals)
        {
            str(s, elem);
            s.append(", ");
        }
        s.pop_back();
        s.back() = ']';
    }
}
template <class A, class B>
void str(std::string& s, std::pair<A, B> const& pair)
{
    s += '(';
    str(s, pair.first);
    s += ", ";
    str(s, pair.second);
    s += ')';
}
class BoundingBox;
void str(std::string& s, BoundingBox const& bb) { s += "BoundingBox"; }

int getMPIRank()
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    return rank;
}

FILE*         get_file() {
    static std::string filename = "rank" + std::to_string(getMPIRank()) + ".log"; 
    static int opened = puts("Opening file"); 
    static FILE* file = fopen(filename.data(), "w");  

    return file;
}
template <class... T>
void mpi_err(T&&... things)
{
    auto file = get_file(); 
    std::string        s;
    char               _[]{(str(s, things), '\0')...};
    s.push_back('\n');
    fwrite(s.data(), 1, s.size(), file);
    fflush(file);
}
