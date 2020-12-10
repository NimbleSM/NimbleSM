target_sources(nimble PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}/parallel_contact_manager.cc
    )

if (HAVE_BVH)
    target_sources(nimble PRIVATE ${CMAKE_CURRENT_LIST_DIR}/bvh_contact_manager.cc)
endif()

if (HAVE_ARBORX AND USE_PURE_MPI)
    set(NIMBLE_PUBLIC_HEADERS ${NIMBLE_PUBLIC_HEADERS}
            ${CMAKE_CURRENT_LIST_DIR}/arborx_parallel_contact_manager.h
            )
    target_sources(nimble PRIVATE
            ${CMAKE_CURRENT_LIST_DIR}/arborx_parallel_contact_manager.cc
        )
endif()
