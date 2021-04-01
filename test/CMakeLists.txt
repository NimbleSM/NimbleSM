#
# Define common variables
#

set(exe_loc ${CMAKE_BINARY_DIR}/src/NimbleSM_MPI)

set(nimble_exe ${exe_loc})
set(nimble_exe_cli_arg "")
set(nimble_exe_flag "")

set(nimble_kokkos ${exe_loc})
set(nimble_kokkos_cli_arg "use_kokkos")
set(nimble_kokkos_flag "-kokkos")

set(nimble_tpetra ${exe_loc})
set(nimble_tpetra_cli_arg "use_tpetra")
set(nimble_tpetra_flag "-tpetra")

#
# Include test directories
#

add_subdirectory(serial)

add_subdirectory(explicit)

if (NIMBLE_HAVE_TRILINOS)
  add_subdirectory(quasistatic)
endif()

if (NIMBLE_HAVE_ARBORX)
  add_subdirectory(contact)
endif()

#
# Copy files
#

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/run_test.py
               ${CMAKE_CURRENT_BINARY_DIR}/run_test.py COPYONLY)

