target_sources(nimble PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}/serial_contact_manager.cc
    )

if (HAVE_ARBORX)
    target_sources(nimble PRIVATE ${CMAKE_CURRENT_LIST_DIR}/arborx_serial_contact_manager.cc)
endif()
