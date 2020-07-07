target_sources(nimble PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}/parallel_contact_manager.cc
    )

if (HAVE_BVH)
    target_sources(nimble PRIVATE ${CMAKE_CURRENT_LIST_DIR}/bvh_contact_manager.cc)
endif()
