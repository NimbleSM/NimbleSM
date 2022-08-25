set(prefix "cube")

foreach (ext "in" "g" "gold.e" "exodiff")
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${prefix}.${ext}
                 ${CMAKE_CURRENT_BINARY_DIR}/${prefix}.${ext} COPYONLY)
endforeach()

add_test(NAME ${cube}-serial
         COMMAND Python::Interpreter ../../run_exodiff_test.py --executable "${nimble_exe}" --cli-flag "" --input-deck ${cube}.in --num-ranks 1
        )

if(NIMBLE_HAVE_KOKKOS)
  add_test(NAME ${cube}-serial-kokkos
           COMMAND Python::Interpreter ../../run_exodiff_test.py --executable "${nimble_exe}" --cli-flag "use_kokkos" --input-deck ${cube}.in --num-ranks 1
          )
endif()
