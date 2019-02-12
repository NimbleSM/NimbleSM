cmake_minimum_required(VERSION 3.1)
project(mpi_buckets)

find_package(MPI REQUIRED)

include_directories(${MPI_INCLUDE_PATH} include)

file(GLOB TESTS "tests/*.cc")

add_executable(run_tests ${TESTS})
target_link_libraries(run_tests ${MPI_LIBRARIES})

