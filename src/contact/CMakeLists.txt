target_sources(nimble PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}/kokkos_contact_manager.cc
    )

if (NIMBLE_HAVE_MPI)
    add_subdirectory(parallel)
endif()