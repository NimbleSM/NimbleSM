#
# Define common variables
#

set(exe_loc ${CMAKE_BINARY_DIR}/src/nimble.exe)

set(nimble_exe ${exe_loc})
set(nimble_exe_cli_arg "")
set(nimble_exe_flag "")

set(nimble_kokkos_cli_arg "use_kokkos")
set(nimble_kokkos_flag "-kokkos")

set(nimble_tpetra_cli_arg "use_tpetra")
set(nimble_tpetra_flag "-tpetra")

#
# Include test directories
#

#--- dynamics
add_subdirectory(brick_with_fibers)
add_subdirectory(notched_plate_native_hypoelastic)
add_subdirectory(notched_plate_native_neohookean)
add_subdirectory(rigid_body_motion)
add_subdirectory(simple_deformation_modes)
add_subdirectory(single_elem_complex_displacement)
add_subdirectory(wave_in_bar)

if (NIMBLE_HAVE_TRILINOS)
# This will become
#    add_subdirectory(quasistatic)
    add_subdirectory(brick_with_fiber_periodic)
endif()

if (NIMBLE_HAVE_KOKKOS)
    add_subdirectory(contact_entity_creation)
endif()

if (NIMBLE_HAVE_ARBORX)
    add_subdirectory(cubes_contact)
    add_subdirectory(sphere_plate_contact)
endif()

#
# Copy files
#

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/run_exodiff_test.py
        ${CMAKE_CURRENT_BINARY_DIR}/run_exodiff_test.py COPYONLY)
