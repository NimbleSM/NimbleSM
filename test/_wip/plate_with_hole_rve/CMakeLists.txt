
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/plate_with_hole_rve.in
               ${CMAKE_CURRENT_BINARY_DIR}/plate_with_hole_rve.in COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/plate_with_hole_rve.g
               ${CMAKE_CURRENT_BINARY_DIR}/plate_with_hole_rve.g COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/plate_with_hole_rve.g.3.0
               ${CMAKE_CURRENT_BINARY_DIR}/plate_with_hole_rve.g.3.0 COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/plate_with_hole_rve.g.3.1
               ${CMAKE_CURRENT_BINARY_DIR}/plate_with_hole_rve.g.3.1 COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/plate_with_hole_rve.g.3.2
               ${CMAKE_CURRENT_BINARY_DIR}/plate_with_hole_rve.g.3.2 COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/rve.g
               ${CMAKE_CURRENT_BINARY_DIR}/rve.g COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/plate_with_hole_rve.gold.e
               ${CMAKE_CURRENT_BINARY_DIR}/plate_with_hole_rve.gold.e COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/plate_with_hole_rve.exodiff
               ${CMAKE_CURRENT_BINARY_DIR}/plate_with_hole_rve.exodiff COPYONLY)

if(NOT NIMBLE_HAVE_KOKKOS)
  add_test(plate_with_hole_rve_serial python ../run_exodiff_test.py --executable "${NimbleSM_Serial}" --input-deck plate_with_hole_rve.in)
endif()

if(NIMBLE_HAVE_DARMA)
  add_test(plate_with_hole_rve_darma_np1 python ../run_exodiff_test.py --executable "${NimbleSM_Darma}" --input-deck plate_with_hole_rve.in --num-ranks 1 --num-virtual-ranks 1)
  add_test(plate_with_hole_rve_darma_np3 python ../run_exodiff_test.py --executable "${NimbleSM_Darma}" --input-deck plate_with_hole_rve.in --num-ranks 2 --num-virtual-ranks 3)
endif()

if(NIMBLE_HAVE_MPI AND NOT NIMBLE_HAVE_KOKKOS)
  add_test(plate_with_hole_rve_mpi_np1 python ../run_exodiff_test.py --executable "${NimbleSM_MPI}" --input-deck plate_with_hole_rve.in --num-ranks 1)
  add_test(plate_with_hole_rve_mpi_np3 python ../run_exodiff_test.py --executable "${NimbleSM_MPI}" --input-deck plate_with_hole_rve.in --num-ranks 3)
endif()

if(NIMBLE_HAVE_TRILINOS)
  add_test(plate_with_hole_rve_tpetra_np1 python ../run_exodiff_test.py --executable "${NimbleSM_Tpetra}" --input-deck plate_with_hole_rve.in --num-ranks 1)
  add_test(plate_with_hole_rve_tpetra_np3 python ../run_exodiff_test.py --executable "${NimbleSM_Tpetra}" --input-deck plate_with_hole_rve.in --num-ranks 3)
endif()

if(NIMBLE_HAVE_QTHREADS)
  add_test(plate_with_hole_rve_qthreads_np1_ns1_nwps1 python ../run_exodiff_test.py --executable "${NimbleSM_Qthreads}" --input-deck plate_with_hole_rve.in --num-ranks 1 --qthreads-num-shepherds 1 --qthreads-num-workers-per-shepherd 1)
  add_test(plate_with_hole_rve_qthreads_np1_ns1_nwps3 python ../run_exodiff_test.py --executable "${NimbleSM_Qthreads}" --input-deck plate_with_hole_rve.in --num-ranks 1 --qthreads-num-shepherds 1 --qthreads-num-workers-per-shepherd 3)
endif()
