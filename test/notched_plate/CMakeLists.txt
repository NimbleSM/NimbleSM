
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/notched_plate.g
               ${CMAKE_CURRENT_BINARY_DIR}/notched_plate.g COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/notched_plate.g.6.0
               ${CMAKE_CURRENT_BINARY_DIR}/notched_plate.g.6.0 COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/notched_plate.g.6.1
               ${CMAKE_CURRENT_BINARY_DIR}/notched_plate.g.6.1 COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/notched_plate.g.6.2
               ${CMAKE_CURRENT_BINARY_DIR}/notched_plate.g.6.2 COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/notched_plate.g.6.3
               ${CMAKE_CURRENT_BINARY_DIR}/notched_plate.g.6.3 COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/notched_plate.g.6.4
               ${CMAKE_CURRENT_BINARY_DIR}/notched_plate.g.6.4 COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/notched_plate.g.6.5
               ${CMAKE_CURRENT_BINARY_DIR}/notched_plate.g.6.5 COPYONLY)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/notched_plate_native_hypoelastic.in
               ${CMAKE_CURRENT_BINARY_DIR}/notched_plate_native_hypoelastic.in COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/notched_plate_native_hypoelastic.gold.e
               ${CMAKE_CURRENT_BINARY_DIR}/notched_plate_native_hypoelastic.gold.e COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/notched_plate_native_hypoelastic.exodiff
               ${CMAKE_CURRENT_BINARY_DIR}/notched_plate_native_hypoelastic.exodiff COPYONLY)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/notched_plate_lame_hypoelastic.in
               ${CMAKE_CURRENT_BINARY_DIR}/notched_plate_lame_hypoelastic.in COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/notched_plate_lame_hypoelastic.gold.e
               ${CMAKE_CURRENT_BINARY_DIR}/notched_plate_lame_hypoelastic.gold.e COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/notched_plate_lame_hypoelastic.exodiff
               ${CMAKE_CURRENT_BINARY_DIR}/notched_plate_lame_hypoelastic.exodiff COPYONLY)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/notched_plate_ngp_lame_hypoelastic.in
               ${CMAKE_CURRENT_BINARY_DIR}/notched_plate_ngp_lame_hypoelastic.in COPYONLY)
execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
                ${CMAKE_CURRENT_BINARY_DIR}/notched_plate_native_hypoelastic.gold.e
                ${CMAKE_CURRENT_BINARY_DIR}/notched_plate_ngp_lame_hypoelastic.gold.e)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/notched_plate_ngp_lame_hypoelastic.exodiff
               ${CMAKE_CURRENT_BINARY_DIR}/notched_plate_ngp_lame_hypoelastic.exodiff COPYONLY)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/notched_plate_native_neohookean.in
               ${CMAKE_CURRENT_BINARY_DIR}/notched_plate_native_neohookean.in COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/notched_plate_native_neohookean.gold.e
               ${CMAKE_CURRENT_BINARY_DIR}/notched_plate_native_neohookean.gold.e COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/notched_plate_native_neohookean.exodiff
               ${CMAKE_CURRENT_BINARY_DIR}/notched_plate_native_neohookean.exodiff COPYONLY)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/notched_plate_lame_neohookean.in
               ${CMAKE_CURRENT_BINARY_DIR}/notched_plate_lame_neohookean.in COPYONLY)
execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
                ${CMAKE_CURRENT_BINARY_DIR}/notched_plate_native_neohookean.gold.e
                ${CMAKE_CURRENT_BINARY_DIR}/notched_plate_lame_neohookean.gold.e)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/notched_plate_lame_neohookean.exodiff
               ${CMAKE_CURRENT_BINARY_DIR}/notched_plate_lame_neohookean.exodiff COPYONLY)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/notched_plate_ngp_lame_neohookean.in
               ${CMAKE_CURRENT_BINARY_DIR}/notched_plate_ngp_lame_neohookean.in COPYONLY)
execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
                ${CMAKE_CURRENT_BINARY_DIR}/notched_plate_native_neohookean.gold.e
                ${CMAKE_CURRENT_BINARY_DIR}/notched_plate_ngp_lame_neohookean.gold.e)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/notched_plate_ngp_lame_neohookean.exodiff
               ${CMAKE_CURRENT_BINARY_DIR}/notched_plate_ngp_lame_neohookean.exodiff COPYONLY)

add_test(notched_plate_native_hypoelastic_serial python ../run_exodiff_test.py --executable "${NimbleSM_Serial}" --input-deck notched_plate_native_hypoelastic.in)
add_test(notched_plate_native_neohookean_serial python ../run_exodiff_test.py --executable "${NimbleSM_Serial}" --input-deck notched_plate_native_neohookean.in)

if(NIMBLE_HAVE_MPI)
 add_test(notched_plate_native_hypoelastic_mpi_np1 python ../run_exodiff_test.py --executable "${NimbleSM_MPI}" --input-deck notched_plate_native_hypoelastic.in --num-ranks 1)
 add_test(notched_plate_native_hypoelastic_mpi_np6 python ../run_exodiff_test.py --executable "${NimbleSM_MPI}" --input-deck notched_plate_native_hypoelastic.in --num-ranks 6)
 add_test(notched_plate_native_neohookean_mpi_np1 python ../run_exodiff_test.py --executable "${NimbleSM_MPI}" --input-deck notched_plate_native_neohookean.in --num-ranks 1)
 add_test(notched_plate_native_neohookean_mpi_np6 python ../run_exodiff_test.py --executable "${NimbleSM_MPI}" --input-deck notched_plate_native_neohookean.in --num-ranks 6)
endif()

if(NIMBLE_HAVE_KOKKOS)
 add_test(notched_plate_native_hypoelastic_kokkos_np1 python ../run_exodiff_test.py --executable "${NimbleSM_Kokkos}" --input-deck notched_plate_native_hypoelastic.in --num-ranks 1)
 add_test(notched_plate_native_hypoelastic_kokkos_np6 python ../run_exodiff_test.py --executable "${NimbleSM_Kokkos}" --input-deck notched_plate_native_hypoelastic.in --num-ranks 6)
 add_test(notched_plate_native_neohookean_kokkos_np1 python ../run_exodiff_test.py --executable "${NimbleSM_Kokkos}" --input-deck notched_plate_native_neohookean.in --num-ranks 1)
 add_test(notched_plate_native_neohookean_kokkos_np6 python ../run_exodiff_test.py --executable "${NimbleSM_Kokkos}" --input-deck notched_plate_native_neohookean.in --num-ranks 6)
endif()

# DEACTIVATED UNTIL ROTATED/UNROTATED STRESS BUG IS FIXED
#if(NIMBLE_HAVE_EXTRAS AND NIMBLE_HAVE_MPI)
# add_test(notched_plate_lame_hypoelastic_mpi_np1 python ../run_exodiff_test.py --executable "${NimbleSM_MPI}" --input-deck notched_plate_lame_hypoelastic.in --num-ranks 1)
# add_test(notched_plate_lame_hypoelastic_mpi_np6 python ../run_exodiff_test.py --executable "${NimbleSM_MPI}" --input-deck notched_plate_lame_hypoelastic.in --num-ranks 6)
# add_test(notched_plate_lame_neohookean_mpi_np1 python ../run_exodiff_test.py --executable "${NimbleSM_MPI}" --input-deck notched_plate_lame_neohookean.in --num-ranks 1)
# add_test(notched_plate_lame_neohookean_mpi_np6 python ../run_exodiff_test.py --executable "${NimbleSM_MPI}" --input-deck notched_plate_lame_neohookean.in --num-ranks 6)
## add_test(notched_plate_ngp_lame_hypoelastic_kokkos_np1 python ../run_exodiff_test.py --executable "${NimbleSM_Kokkos}" --input-deck notched_plate_ngp_lame_hypoelastic.in --num-ranks 1)
## add_test(notched_plate_ngp_lame_hypoelastic_kokkos_np6 python ../run_exodiff_test.py --executable "${NimbleSM_Kokkos}" --input-deck notched_plate_ngp_lame_hypoelastic.in --num-ranks 6)
# add_test(notched_plate_ngp_lame_neohookean_kokkos_np1 python ../run_exodiff_test.py --executable "${NimbleSM_Kokkos}" --input-deck notched_plate_ngp_lame_neohookean.in --num-ranks 1)
# add_test(notched_plate_ngp_lame_neohookean_kokkos_np6 python ../run_exodiff_test.py --executable "${NimbleSM_Kokkos}" --input-deck notched_plate_ngp_lame_neohookean.in --num-ranks 6)
#endif()

if(NIMBLE_HAVE_TRILINOS)
  add_test(notched_plate_native_hypoelastic_tpetra_np1 python ../run_exodiff_test.py --executable "${NimbleSM_Tpetra}" --input-deck notched_plate_native_hypoelastic.in --num-ranks 1)
  add_test(notched_plate_native_hypoelastic_tpetra_np6 python ../run_exodiff_test.py --executable "${NimbleSM_Tpetra}" --input-deck notched_plate_native_hypoelastic.in --num-ranks 6)
  add_test(notched_plate_native_neohookean_tpetra_np1 python ../run_exodiff_test.py --executable "${NimbleSM_Tpetra}" --input-deck notched_plate_native_neohookean.in --num-ranks 1)
  add_test(notched_plate_native_neohookean_tpetra_np6 python ../run_exodiff_test.py --executable "${NimbleSM_Tpetra}" --input-deck notched_plate_native_neohookean.in --num-ranks 6)
endif()

if(NIMBLE_HAVE_QTHREADS)
  add_test(notched_plate_native_hypoelastic_qthreads_np1_ns1_nwps1 python ../run_exodiff_test.py --executable "${NimbleSM_Qthreads}" --input-deck notched_plate_native_hypoelastic.in --num-ranks 1 --qthreads-num-shepherds 1 --qthreads-num-workers-per-shepherd 1)
  add_test(notched_plate_native_hypoelastic_qthreads_np2_ns3_nwps1 python ../run_exodiff_test.py --executable "${NimbleSM_Qthreads}" --input-deck notched_plate_native_hypoelastic.in --num-ranks 2 --qthreads-num-shepherds 3 --qthreads-num-workers-per-shepherd 1)
  add_test(notched_plate_native_neohookean_qthreads_np1_ns1_nwps1 python ../run_exodiff_test.py --executable "${NimbleSM_Qthreads}" --input-deck notched_plate_native_neohookean.in --num-ranks 1 --qthreads-num-shepherds 1 --qthreads-num-workers-per-shepherd 1)
  add_test(notched_plate_native_neohookean_qthreads_np2_ns3_nwps1 python ../run_exodiff_test.py --executable "${NimbleSM_Qthreads}" --input-deck notched_plate_native_neohookean.in --num-ranks 2 --qthreads-num-shepherds 3 --qthreads-num-workers-per-shepherd 1)
endif()

if(NIMBLE_HAVE_DARMA)
 add_test(notched_plate_native_hypoelastic_darma_np1 python ../run_exodiff_test.py --executable "${NimbleSM_Darma}" --input-deck notched_plate_native_hypoelastic.in --num-ranks 1 --num-virtual-ranks 1 --charm-run ${charmrun})
 add_test(notched_plate_native_hypoelastic_darma_np6 python ../run_exodiff_test.py --executable "${NimbleSM_Darma}" --input-deck notched_plate_native_hypoelastic.in --num-ranks 2 --num-virtual-ranks 6 --charm-run ${charmrun})
 add_test(notched_plate_native_neohookean_darma_np1 python ../run_exodiff_test.py --executable "${NimbleSM_Darma}" --input-deck notched_plate_native_neohookean.in --num-ranks 1 --num-virtual-ranks 1 --charm-run ${charmrun})
 add_test(notched_plate_native_neohookean_darma_np6 python ../run_exodiff_test.py --executable "${NimbleSM_Darma}" --input-deck notched_plate_native_neohookean.in --num-ranks 2 --num-virtual-ranks 6 --charm-run ${charmrun})
endif()
