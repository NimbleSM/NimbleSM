#
# Define common variables
#

set(exe_loc ${CMAKE_BINARY_DIR}/src/NimbleSM)

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

#--- contact
add_subdirectory(contact)

#--- dynamics
add_subdirectory(dynamics)

#--- quasistatic
add_subdirectory(quasistatic)

#
# Copy files
#

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/run_exodiff_test.py
        ${CMAKE_CURRENT_BINARY_DIR}/run_exodiff_test.py COPYONLY)
