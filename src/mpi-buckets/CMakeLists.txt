cmake_minimum_required(VERSION 3.1)
project(mpi_buckets)

set(CMAKE_CXX_FLAGS "-std=c++11 -O3")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
find_package(MPI REQUIRED)

include_directories(${MPI_INCLUDE_PATH} src)

file(GLOB TESTS "tests/*.cc")

add_executable(run_tests ${TESTS})
target_link_libraries(run_tests ${MPI_LIBRARIES})
