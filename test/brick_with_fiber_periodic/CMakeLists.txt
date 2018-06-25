
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/brick_with_fiber_periodic.in
               ${CMAKE_CURRENT_BINARY_DIR}/brick_with_fiber_periodic.in COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/brick_with_fiber_periodic.g
               ${CMAKE_CURRENT_BINARY_DIR}/brick_with_fiber_periodic.g COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/brick_with_fiber_periodic.gold.e
               ${CMAKE_CURRENT_BINARY_DIR}/brick_with_fiber_periodic.gold.e COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/brick_with_fiber_periodic.exodiff
               ${CMAKE_CURRENT_BINARY_DIR}/brick_with_fiber_periodic.exodiff COPYONLY)

add_test(brick_with_fiber_periodic_serial python ../run_exodiff_test.py --executable "${NimbleSM_Serial}" --input-deck brick_with_fiber_periodic.in)

if(NIMBLE_HAVE_TRILINOS)
  add_test(brick_with_fiber_periodic_tpetra_np1 python ../run_exodiff_test.py --executable "${NimbleSM_Tpetra}" --input-deck brick_with_fiber_periodic.in --num-ranks 1)
endif()
