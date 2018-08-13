
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/brick_with_fibers.in
               ${CMAKE_CURRENT_BINARY_DIR}/brick_with_fibers.in COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/brick_with_fibers.g
               ${CMAKE_CURRENT_BINARY_DIR}/brick_with_fibers.g COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/brick_with_fibers.g.2.0
               ${CMAKE_CURRENT_BINARY_DIR}/brick_with_fibers.g.2.0 COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/brick_with_fibers.g.2.1
               ${CMAKE_CURRENT_BINARY_DIR}/brick_with_fibers.g.2.1 COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/brick_with_fibers.gold.e
               ${CMAKE_CURRENT_BINARY_DIR}/brick_with_fibers.gold.e COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/brick_with_fibers.exodiff
               ${CMAKE_CURRENT_BINARY_DIR}/brick_with_fibers.exodiff COPYONLY)

add_test(brick_with_fibers_serial python ../run_exodiff_test.py --executable "${NimbleSM_Serial}" --input-deck brick_with_fibers.in)

if(NIMBLE_HAVE_DARMA)
  add_test(brick_with_fibers_darma_np1 python ../run_exodiff_test.py --executable "${NimbleSM_Darma}" --input-deck brick_with_fibers.in --num-ranks 1 --num-virtual-ranks 1 --charm-run ${charmrun})
  add_test(brick_with_fibers_darma_np2 python ../run_exodiff_test.py --executable "${NimbleSM_Darma}" --input-deck brick_with_fibers.in --num-ranks 2 --num-virtual-ranks 2 --charm-run ${charmrun})
endif()

if(NIMBLE_HAVE_MPI)
  add_test(brick_with_fibers_mpi_np1 python ../run_exodiff_test.py --executable "${NimbleSM_MPI}" --input-deck brick_with_fibers.in --num-ranks 1)
  add_test(brick_with_fibers_mpi_np2 python ../run_exodiff_test.py --executable "${NimbleSM_MPI}" --input-deck brick_with_fibers.in --num-ranks 2)
endif()

if(NIMBLE_HAVE_KOKKOS)
  add_test(brick_with_fibers_kokkos_np1 python ../run_exodiff_test.py --executable "${NimbleSM_Kokkos}" --input-deck brick_with_fibers.in --num-ranks 1)
  add_test(brick_with_fibers_kokkos_np2 python ../run_exodiff_test.py --executable "${NimbleSM_Kokkos}" --input-deck brick_with_fibers.in --num-ranks 2)
endif()

if(NIMBLE_HAVE_TRILINOS)
  add_test(brick_with_fibers_tpetra_np1 python ../run_exodiff_test.py --executable "${NimbleSM_Tpetra}" --input-deck brick_with_fibers.in --num-ranks 1)
  add_test(brick_with_fibers_tpetra_np2 python ../run_exodiff_test.py --executable "${NimbleSM_Tpetra}" --input-deck brick_with_fibers.in --num-ranks 2)
endif()

if(NIMBLE_HAVE_QTHREADS)
  add_test(brick_with_fibers_qthreads_np1_ns1_nwps1 python ../run_exodiff_test.py --executable "${NimbleSM_Qthreads}" --input-deck brick_with_fibers.in --num-ranks 1 --qthreads-num-shepherds 1 --qthreads-num-workers-per-shepherd 1)
  add_test(brick_with_fibers_qthreads_np1_ns1_nwps2 python ../run_exodiff_test.py --executable "${NimbleSM_Qthreads}" --input-deck brick_with_fibers.in --num-ranks 1 --qthreads-num-shepherds 1 --qthreads-num-workers-per-shepherd 2)
endif()
