
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/single_elem_complex_displacement.in
               ${CMAKE_CURRENT_BINARY_DIR}/single_elem_complex_displacement.in COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/single_elem_complex_displacement.g
               ${CMAKE_CURRENT_BINARY_DIR}/single_elem_complex_displacement.g COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/single_elem_complex_displacement.gold.e
               ${CMAKE_CURRENT_BINARY_DIR}/single_elem_complex_displacement.gold.e COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/single_elem_complex_displacement.exodiff
               ${CMAKE_CURRENT_BINARY_DIR}/single_elem_complex_displacement.exodiff COPYONLY)

add_test(single_elem_complex_displacement_serial python ../run_exodiff_test.py --executable "${NimbleSM_Serial}" --input-deck single_elem_complex_displacement.in)

if(NIMBLE_HAVE_MPI)
  add_test(single_elem_complex_displacement_mpi_np1 python ../run_exodiff_test.py --executable "${NimbleSM_MPI}" --input-deck single_elem_complex_displacement.in --num-ranks 1)
endif()

if(NIMBLE_HAVE_KOKKOS)
  add_test(single_elem_complex_displacement_kokkos_np1 python ../run_exodiff_test.py --executable "${NimbleSM_Kokkos}" --input-deck single_elem_complex_displacement.in --num-ranks 1)
endif()

if(NIMBLE_HAVE_TRILINOS)
  add_test(single_elem_complex_displacement_tpetra_np1 python ../run_exodiff_test.py --executable "${NimbleSM_Tpetra}" --input-deck single_elem_complex_displacement.in --num-ranks 1)
endif()

if(NIMBLE_HAVE_QTHREADS)
  add_test(single_elem_complex_displacement_qthreads_np1_ns1_nwps1 python ../run_exodiff_test.py --executable "${NimbleSM_Qthreads}" --input-deck single_elem_complex_displacement.in --num-ranks 1 --qthreads-num-shepherds 1 --qthreads-num-workers-per-shepherd 1)
endif()
