add_executable(preskin main.cpp boundary.cpp)
target_include_directories(preskin PRIVATE
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/tpl>
        ${SEACASIoss_INCLUDE_DIRS})
target_link_libraries(preskin PRIVATE Kokkos::kokkos)
target_link_libraries(preskin PRIVATE SEACASIoss::Ioss SEACASIoss::Ionit)
target_link_Libraries(preskin PRIVATE MPI::MPI_CXX)