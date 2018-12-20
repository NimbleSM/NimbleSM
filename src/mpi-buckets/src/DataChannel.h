#pragma once 
#include <mpi.h>

struct ActionGroup {
    constexpr static auto any_source = MPI_ANY_SOURCE;

    MPI_Comm              comm;
    int                   tag;

    MPI_Request Irecv_any(char* buff, size_t count) {
        MPI_Request request;
        MPI_Irecv((void*)buff, count, MPI_BYTE, any_source, tag, comm, &request); 
        return request;
    }
};
