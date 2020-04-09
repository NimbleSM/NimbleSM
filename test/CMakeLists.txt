
if(NIMBLE_HAVE_CHARM)
  set(charmrun ${PROJECT_BINARY_DIR}/src/charmrun)
else()
  set(charmrun "not available")
endif()

set(NimbleSM_Serial ${PROJECT_BINARY_DIR}/src/NimbleSM_Serial)
set(NimbleSM_MPI ${PROJECT_BINARY_DIR}/src/NimbleSM_MPI)
set(NimbleSM_Kokkos ${PROJECT_BINARY_DIR}/src/NimbleSM_Kokkos)
set(NimbleSM_Darma ${PROJECT_BINARY_DIR}/src/NimbleSM_Darma)
set(NimbleSM_Qthreads ${PROJECT_BINARY_DIR}/src/NimbleSM_Qthreads)
set(NimbleSM_Tpetra ${PROJECT_BINARY_DIR}/src/NimbleSM_Tpetra)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/run_exodiff_test.py
               ${CMAKE_CURRENT_BINARY_DIR}/run_exodiff_test.py COPYONLY)

add_subdirectory(rigid_body_motion)
add_subdirectory(simple_deformation_modes)
add_subdirectory(wave_in_bar)
add_subdirectory(notched_plate)
add_subdirectory(brick_with_fibers)
add_subdirectory(brick_with_fiber_periodic)
add_subdirectory(quasistatic_tension)
add_subdirectory(quasistatic_tension_J2_lame)
add_subdirectory(contact_entity_creation)
add_subdirectory(two_body_contact)

if(NIMBLE_HAVE_UQ)
  add_subdirectory(uq_wave_in_bar)
endif()

# NOTE failing
#add_subdirectory(cylinder_with_hole_multiblock_rve)
#add_subdirectory(plate_with_hole_rve)
#add_subdirectory(tiny_fe_squared)
