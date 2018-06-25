
if(NIMBLE_HAVE_CHARM)
  set(charmrun ${PROJECT_BINARY_DIR}/src/charmrun)
else()
  set(charmrun "not available")
endif()

set(NimbleSM_Serial ${PROJECT_BINARY_DIR}/src/NimbleSM_Serial)
set(NimbleSM_Tpetra ${PROJECT_BINARY_DIR}/src/NimbleSM_Tpetra)
set(NimbleSM_Darma ${PROJECT_BINARY_DIR}/src/NimbleSM_Darma)
set(NimbleSM_MPI ${PROJECT_BINARY_DIR}/src/NimbleSM_MPI)
set(NimbleSM_Qthreads ${PROJECT_BINARY_DIR}/src/NimbleSM_Qthreads)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/run_exodiff_test.py
               ${CMAKE_CURRENT_BINARY_DIR}/run_exodiff_test.py COPYONLY)

add_subdirectory(rigid_body_motion)
add_subdirectory(simple_deformation_modes)
add_subdirectory(wave_in_bar)
add_subdirectory(notched_plate)
add_subdirectory(notched_plate_lame)
add_subdirectory(notched_plate_hypoelastic_lame)
add_subdirectory(cylinder_with_hole_multiblock_rve)
add_subdirectory(plate_with_hole_rve)
add_subdirectory(brick_with_fibers)
add_subdirectory(brick_with_fiber_periodic)
add_subdirectory(quasistatic_tension)
add_subdirectory(quasistatic_tension_J2_lame)
add_subdirectory(tiny_fe_squared)
