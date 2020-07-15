
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/quasistatic_tension.in
               ${CMAKE_CURRENT_BINARY_DIR}/quasistatic_tension.in COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/quasistatic_tension.g
               ${CMAKE_CURRENT_BINARY_DIR}/quasistatic_tension.g COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/quasistatic_tension.gold.e
               ${CMAKE_CURRENT_BINARY_DIR}/quasistatic_tension.gold.e COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/quasistatic_tension.exodiff
               ${CMAKE_CURRENT_BINARY_DIR}/quasistatic_tension.exodiff COPYONLY)

add_test(quasistatic_tension_serial python ../run_exodiff_test.py  --executable "${NimbleSM_Serial}" --input-deck quasistatic_tension.in)

#
# Test with Trilinos (=TPetra) fails
#
#if(NIMBLE_HAVE_TRILINOS)
#  add_test(quasistatic_tension_tpetra_np1 python ../run_exodiff_test.py --executable "${NimbleSM_Tpetra}" --input-deck quasistatic_tension.in --num-ranks 1)
#endif()
#
