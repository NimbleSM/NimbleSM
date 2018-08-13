
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/simple_deformation_modes.in
               ${CMAKE_CURRENT_BINARY_DIR}/simple_deformation_modes.in COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/simple_deformation_modes.g
               ${CMAKE_CURRENT_BINARY_DIR}/simple_deformation_modes.g COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/simple_deformation_modes.gold.e
               ${CMAKE_CURRENT_BINARY_DIR}/simple_deformation_modes.gold.e COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/simple_deformation_modes.exodiff
               ${CMAKE_CURRENT_BINARY_DIR}/simple_deformation_modes.exodiff COPYONLY)

add_test(simple_deformation_modes_serial python ../run_exodiff_test.py --executable "${NimbleSM_Serial}" --input-deck simple_deformation_modes.in)

if(NIMBLE_HAVE_DARMA)
  add_test(simple_deformation_modes_darma_np1 python ../run_exodiff_test.py --executable "${NimbleSM_Darma}" --input-deck simple_deformation_modes.in --num-ranks 1 --num-virtual-ranks 1 --charm-run ${charmrun})
endif()

if(NIMBLE_HAVE_MPI)
  add_test(simple_deformation_modes_mpi_np1 python ../run_exodiff_test.py --executable "${NimbleSM_MPI}" --input-deck simple_deformation_modes.in --num-ranks 1)
endif()

if(NIMBLE_HAVE_KOKKOS)
  add_test(simple_deformation_modes_kokkos_np1 python ../run_exodiff_test.py --executable "${NimbleSM_Kokkos}" --input-deck simple_deformation_modes.in --num-ranks 1)
endif()

if(NIMBLE_HAVE_TRILINOS)
  add_test(simple_deformation_modes_tpetra_np1 python ../run_exodiff_test.py --executable "${NimbleSM_Tpetra}" --input-deck simple_deformation_modes.in --num-ranks 1)
endif()

if(NIMBLE_HAVE_QTHREADS)
  add_test(simple_deformation_modes_qthreads_np1_ns1_nwps1 python ../run_exodiff_test.py --executable "${NimbleSM_Qthreads}" --input-deck simple_deformation_modes.in --num-ranks 1 --qthreads-num-shepherds 1 --qthreads-num-workers-per-shepherd 1)
endif()
