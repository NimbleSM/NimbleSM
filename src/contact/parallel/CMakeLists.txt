target_sources(nimble PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}/parallel_contact_manager.cc
    )

if (NIMBLE_HAVE_BVH)
    target_sources(nimble PRIVATE ${CMAKE_CURRENT_LIST_DIR}/bvh_contact_manager.cc)
endif()

if (NIMBLE_HAVE_ARBORX AND NIMBLE_HAVE_MPI)
    set(NIMBLE_PUBLIC_HEADERS ${NIMBLE_PUBLIC_HEADERS}
            ${CMAKE_CURRENT_LIST_DIR}/arborx_parallel_contact_manager.h
            )
    target_sources(nimble PRIVATE
            ${CMAKE_CURRENT_LIST_DIR}/arborx_parallel_contact_manager.cc
        )
endif()
