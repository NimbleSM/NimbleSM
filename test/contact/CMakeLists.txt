#
# This series of tests need a contact manager library
# (for example, ArborX)
#
if (NIMBLE_HAVE_ARBORX)
  add_subdirectory(cubes_contact)
  add_subdirectory(sphere_plate_contact)
endif()

