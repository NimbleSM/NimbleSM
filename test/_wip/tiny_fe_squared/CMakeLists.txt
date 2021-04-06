
configure_file(tiny_fe_squared.in
               ${CMAKE_CURRENT_BINARY_DIR}/tiny_fe_squared.in COPYONLY)
configure_file(tiny_fe_squared.g
               ${CMAKE_CURRENT_BINARY_DIR}/tiny_fe_squared.g COPYONLY)
configure_file(tiny_fe_squared.g.2.0
               ${CMAKE_CURRENT_BINARY_DIR}/tiny_fe_squared.g.2.0 COPYONLY)
configure_file(tiny_fe_squared.g.2.1
               ${CMAKE_CURRENT_BINARY_DIR}/tiny_fe_squared.g.2.1 COPYONLY)
configure_file(tiny_fe_squared_rve.g
               ${CMAKE_CURRENT_BINARY_DIR}/tiny_fe_squared_rve.g COPYONLY)
configure_file(tiny_fe_squared.gold.e
               ${CMAKE_CURRENT_BINARY_DIR}/tiny_fe_squared.gold.e COPYONLY)
configure_file(tiny_fe_squared.exodiff
               ${CMAKE_CURRENT_BINARY_DIR}/tiny_fe_squared.exodiff COPYONLY)

IF(NOT NIMBLE_HAVE_KOKKOS)
  add_test(tiny_fe_squared_serial python ../../run_exodiff_test.py --executable "${NimbleSM_Serial}" --input-deck tiny_fe_squared.in)
ENDIF()

if(NIMBLE_HAVE_DARMA)
  add_test(tiny_fe_squared_darma_np1 python ../../run_exodiff_test.py --executable "${NimbleSM_Darma}" --input-deck tiny_fe_squared.in --num-ranks 1 --num-virtual-ranks 1)
  add_test(tiny_fe_squared_darma_np2 python ../../run_exodiff_test.py --executable "${NimbleSM_Darma}" --input-deck tiny_fe_squared.in --num-ranks 2 --num-virtual-ranks 2)
endif()

if(NIMBLE_HAVE_MPI AND NOT NIMBLE_HAVE_KOKKOS)
  add_test(tiny_fe_squared_mpi_np1 python ../../run_exodiff_test.py --executable "${NimbleSM_MPI}" --input-deck tiny_fe_squared.in --num-ranks 1)
  add_test(tiny_fe_squared_mpi_np2 python ../../run_exodiff_test.py --executable "${NimbleSM_MPI}" --input-deck tiny_fe_squared.in --num-ranks 2)
endif()

if(NIMBLE_HAVE_TRILINOS)
  add_test(tiny_fe_squared_tpetra_np1 python ../../run_exodiff_test.py --executable "${NimbleSM_Tpetra}" --input-deck tiny_fe_squared.in --num-ranks 1)
  add_test(tiny_fe_squared_tpetra_np2 python ../../run_exodiff_test.py --executable "${NimbleSM_Tpetra}" --input-deck tiny_fe_squared.in --num-ranks 2)
endif()

if(NIMBLE_HAVE_QTHREADS)
  add_test(tiny_fe_squared_qthreads_np1_ns1_nwps1 python ../../run_exodiff_test.py --executable "${NimbleSM_Qthreads}" --input-deck tiny_fe_squared.in --num-ranks 1 --qthreads-num-shepherds 1 --qthreads-num-workers-per-shepherd 1)
  add_test(tiny_fe_squared_qthreads_np1_ns2_nwps1 python ../../run_exodiff_test.py --executable "${NimbleSM_Qthreads}" --input-deck tiny_fe_squared.in --num-ranks 1 --qthreads-num-shepherds 2 --qthreads-num-workers-per-shepherd 1)
endif()
