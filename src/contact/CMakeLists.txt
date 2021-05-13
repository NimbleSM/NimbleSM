
if (NIMBLE_HAVE_KOKKOS)
    set(NIMBLE_PUBLIC_HEADERS ${NIMBLE_PUBLIC_HEADERS}
            ${CMAKE_CURRENT_LIST_DIR}/kokkos_contact_manager.h
        )
    target_sources(nimble PRIVATE
            ${CMAKE_CURRENT_LIST_DIR}/kokkos_contact_manager.cc
        )
endif()

if (NIMBLE_HAVE_MPI)
    add_subdirectory(parallel)
endif()

add_subdirectory(serial)

if (NIMBLE_HAVE_ARBORX)
    set(NIMBLE_PUBLIC_HEADERS ${NIMBLE_PUBLIC_HEADERS}
            ${CMAKE_CURRENT_LIST_DIR}/arborx_utils.h
            )
endif()

