
set(NIMBLE_PUBLIC_HEADERS ${NIMBLE_PUBLIC_HEADERS}
    ${CMAKE_CURRENT_LIST_DIR}/serial_contact_manager.h
    )

target_sources(nimble PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}/serial_contact_manager.cc
    )

if (HAVE_ARBORX)
    set(NIMBLE_PUBLIC_HEADERS ${NIMBLE_PUBLIC_HEADERS}
            ${CMAKE_CURRENT_LIST_DIR}/arborx_serial_contact_manager.h
            )
    target_sources(nimble PRIVATE
            ${CMAKE_CURRENT_LIST_DIR}/arborx_serial_contact_manager.cc
        )
endif()
