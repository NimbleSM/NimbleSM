
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/rigid_body_motion.in
               ${CMAKE_CURRENT_BINARY_DIR}/rigid_body_motion.in COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/rigid_body_motion.g
               ${CMAKE_CURRENT_BINARY_DIR}/rigid_body_motion.g COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/rigid_body_motion.g.2.0
               ${CMAKE_CURRENT_BINARY_DIR}/rigid_body_motion.g.2.0 COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/rigid_body_motion.g.2.1
               ${CMAKE_CURRENT_BINARY_DIR}/rigid_body_motion.g.2.1 COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/rigid_body_motion.gold.e
               ${CMAKE_CURRENT_BINARY_DIR}/rigid_body_motion.gold.e COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/rigid_body_motion.exodiff
               ${CMAKE_CURRENT_BINARY_DIR}/rigid_body_motion.exodiff COPYONLY)

add_test(rigid_body_motion_serial python ../run_exodiff_test.py --executable "${NimbleSM_Serial}" --input-deck rigid_body_motion.in)

if(NIMBLE_HAVE_DARMA)
  add_test(rigid_body_motion_darma_np1 python ../run_exodiff_test.py --executable "${NimbleSM_Darma}" --input-deck rigid_body_motion.in --num-ranks 1 --num-virtual-ranks 1 --charm-run ${charmrun})
  add_test(rigid_body_motion_darma_np2 python ../run_exodiff_test.py --executable "${NimbleSM_Darma}" --input-deck rigid_body_motion.in --num-ranks 2 --num-virtual-ranks 2 --charm-run ${charmrun})
endif()

if(NIMBLE_HAVE_MPI)
  add_test(rigid_body_motion_mpi_np1 python ../run_exodiff_test.py --executable "${NimbleSM_MPI}" --input-deck rigid_body_motion.in --num-ranks 1)
  add_test(rigid_body_motion_mpi_np2 python ../run_exodiff_test.py --executable "${NimbleSM_MPI}" --input-deck rigid_body_motion.in --num-ranks 2)
endif()

if(NIMBLE_HAVE_ARBORX)
  add_test(rigid_body_motion_arborx_np1 python ../run_exodiff_test.py --executable "${NimbleSM_ArborX}" --input-deck rigid_body_motion.in --num-ranks 1)
endif()

if(NIMBLE_HAVE_KOKKOS AND NOT NIMBLE_HAVE_ARBORX)
  add_test(rigid_body_motion_kokkos_np1 python ../run_exodiff_test.py --executable "${NimbleSM_Kokkos}" --input-deck rigid_body_motion.in --num-ranks 1)
  if (NIMBLE_HAVE_MPI)
    add_test(rigid_body_motion_kokkos_np2 python ../run_exodiff_test.py --executable "${NimbleSM_Kokkos}" --input-deck rigid_body_motion.in --num-ranks 2)
  endif()
endif()

if(NIMBLE_HAVE_TRILINOS)
  add_test(rigid_body_motion_tpetra_np1 python ../run_exodiff_test.py --executable "${NimbleSM_Tpetra}" --input-deck rigid_body_motion.in --num-ranks 1)
  add_test(rigid_body_motion_tpetra_np2 python ../run_exodiff_test.py --executable "${NimbleSM_Tpetra}" --input-deck rigid_body_motion.in --num-ranks 2)
endif()

if(NIMBLE_HAVE_QTHREADS)
  add_test(rigid_body_motion_qthreads_np1_ns1_nwps1 python ../run_exodiff_test.py --executable "${NimbleSM_Qthreads}" --input-deck rigid_body_motion.in --num-ranks 1 --qthreads-num-shepherds 1 --qthreads-num-workers-per-shepherd 1)
  add_test(rigid_body_motion_qthreads_np1_ns1_nwps2 python ../run_exodiff_test.py --executable "${NimbleSM_Qthreads}" --input-deck rigid_body_motion.in --num-ranks 1 --qthreads-num-shepherds 1 --qthreads-num-workers-per-shepherd 2)
endif()
