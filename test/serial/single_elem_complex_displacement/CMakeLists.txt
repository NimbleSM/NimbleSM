
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/single_elem_complex_displacement.in
               ${CMAKE_CURRENT_BINARY_DIR}/single_elem_complex_displacement.in COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/single_elem_complex_displacement.g
               ${CMAKE_CURRENT_BINARY_DIR}/single_elem_complex_displacement.g COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/single_elem_complex_displacement.gold.e
               ${CMAKE_CURRENT_BINARY_DIR}/single_elem_complex_displacement.gold.e COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/single_elem_complex_displacement.exodiff
               ${CMAKE_CURRENT_BINARY_DIR}/single_elem_complex_displacement.exodiff COPYONLY)

add_test(NAME single_elem_complex_displacement-serial
         COMMAND python ../../run_test.py --executable "${nimble_exe}" --cli-flag "" --input-deck single_elem_complex_displacement.in --num-ranks 1
        )

if(NIMBLE_HAVE_KOKKOS)
  add_test(NAME single_elem_complex_displacement-serial-kokkos
           COMMAND python ../../run_test.py --executable "${nimble_kokkos}" --cli-flag "use_kokkos" --input-deck single_elem_complex_displacement.in --num-ranks 1
          )
endif()

