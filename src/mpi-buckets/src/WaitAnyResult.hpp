#pragma once
#include <mpi.h>

struct WaitAnyResult {
    MPI_Status  status;
    MPI_Request request;
};
