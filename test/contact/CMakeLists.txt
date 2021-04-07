#
# Include test directories
#

if (NIMBLE_HAVE_KOKKOS)
    add_subdirectory(contact_entity_creation)
endif()

if (NIMBLE_HAVE_ARBORX)
    add_subdirectory(cubes_contact)
    add_subdirectory(sphere_plate_contact)
endif()

