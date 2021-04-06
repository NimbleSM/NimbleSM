
configure_file(cylinder_with_hole.in
               ${CMAKE_CURRENT_BINARY_DIR}/cylinder_with_hole.in COPYONLY)
configure_file(cylinder_with_hole.g
               ${CMAKE_CURRENT_BINARY_DIR}/cylinder_with_hole.g COPYONLY)
configure_file(cylinder_with_hole.g.4.0
               ${CMAKE_CURRENT_BINARY_DIR}/cylinder_with_hole.g.4.0 COPYONLY)
configure_file(cylinder_with_hole.g.4.1
               ${CMAKE_CURRENT_BINARY_DIR}/cylinder_with_hole.g.4.1 COPYONLY)
configure_file(cylinder_with_hole.g.4.2
               ${CMAKE_CURRENT_BINARY_DIR}/cylinder_with_hole.g.4.2 COPYONLY)
configure_file(cylinder_with_hole.g.4.3
               ${CMAKE_CURRENT_BINARY_DIR}/cylinder_with_hole.g.4.3 COPYONLY)
configure_file(cylinder_with_hole.g.8.0
               ${CMAKE_CURRENT_BINARY_DIR}/cylinder_with_hole.g.8.0 COPYONLY)
configure_file(cylinder_with_hole.g.8.1
               ${CMAKE_CURRENT_BINARY_DIR}/cylinder_with_hole.g.8.1 COPYONLY)
configure_file(cylinder_with_hole.g.8.2
               ${CMAKE_CURRENT_BINARY_DIR}/cylinder_with_hole.g.8.2 COPYONLY)
configure_file(cylinder_with_hole.g.8.3
               ${CMAKE_CURRENT_BINARY_DIR}/cylinder_with_hole.g.8.3 COPYONLY)
configure_file(cylinder_with_hole.g.8.4
               ${CMAKE_CURRENT_BINARY_DIR}/cylinder_with_hole.g.8.4 COPYONLY)
configure_file(cylinder_with_hole.g.8.5
               ${CMAKE_CURRENT_BINARY_DIR}/cylinder_with_hole.g.8.5 COPYONLY)
configure_file(cylinder_with_hole.g.8.6
               ${CMAKE_CURRENT_BINARY_DIR}/cylinder_with_hole.g.8.6 COPYONLY)
configure_file(cylinder_with_hole.g.8.7
               ${CMAKE_CURRENT_BINARY_DIR}/cylinder_with_hole.g.8.7 COPYONLY)
configure_file(rve.g
               ${CMAKE_CURRENT_BINARY_DIR}/rve.g COPYONLY)
configure_file(cylinder_with_hole.gold.e
               ${CMAKE_CURRENT_BINARY_DIR}/cylinder_with_hole.gold.e COPYONLY)
configure_file(cylinder_with_hole.exodiff
               ${CMAKE_CURRENT_BINARY_DIR}/cylinder_with_hole.exodiff COPYONLY)

if(NOT NIMBLE_HAVE_KOKKOS)
  add_test(cylinder_with_hole_multiblock_rve_serial python ../../run_exodiff_test.py --executable "${NimbleSM_Serial}" --input-deck cylinder_with_hole.in)
endif()

if(NIMBLE_HAVE_DARMA)
  add_test(cylinder_with_hole_multiblock_rve_darma_np1 python ../../run_exodiff_test.py --executable "${NimbleSM_Darma}" --input-deck cylinder_with_hole.in --num-ranks 1 --num-virtual-ranks 1)
  add_test(cylinder_with_hole_multiblock_rve_darma_np4 python ../../run_exodiff_test.py --executable "${NimbleSM_Darma}" --input-deck cylinder_with_hole.in --num-ranks 4 --num-virtual-ranks 8)
endif()

if(NIMBLE_HAVE_MPI AND NOT NIMBLE_HAVE_KOKKOS)
  add_test(cylinder_with_hole_multiblock_rve_mpi_np1 python ../../run_exodiff_test.py --executable "${NimbleSM_MPI}" --input-deck cylinder_with_hole.in --num-ranks 1)
  add_test(cylinder_with_hole_multiblock_rve_mpi_np4 python ../../run_exodiff_test.py --executable "${NimbleSM_MPI}" --input-deck cylinder_with_hole.in --num-ranks 4)
endif()

if(NIMBLE_HAVE_TRILINOS)
  add_test(cylinder_with_hole_multiblock_rve_tpetra_np1 python ../../run_exodiff_test.py --executable "${NimbleSM_Tpetra}" --input-deck cylinder_with_hole.in --num-ranks 1)
  add_test(cylinder_with_hole_multiblock_rve_tpetra_np4 python ../../run_exodiff_test.py --executable "${NimbleSM_Tpetra}" --input-deck cylinder_with_hole.in --num-ranks 4)
endif()

if(NIMBLE_HAVE_QTHREADS)
  add_test(cylinder_with_hole_multiblock_rve_qthreads_np1_ns1_nwps1 python ../../run_exodiff_test.py --executable "${NimbleSM_Qthreads}" --input-deck cylinder_with_hole.in --num-ranks 1 --qthreads-num-shepherds 1 --qthreads-num-workers-per-shepherd 1)
  add_test(cylinder_with_hole_multiblock_rve_qthreads_np2_ns2_nwps1 python ../../run_exodiff_test.py --executable "${NimbleSM_Qthreads}" --input-deck cylinder_with_hole.in --num-ranks 2 --qthreads-num-shepherds 2 --qthreads-num-workers-per-shepherd 1)
endif()
