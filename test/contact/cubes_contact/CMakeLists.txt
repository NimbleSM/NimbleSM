
set(prefix "cubes_contact")

if (NIMBLE_HAVE_ARBORX)
  #
  # When ArborX is set, Kokkos is also present
  #
  foreach (ext "in" "g" "gold.e" "exodiff")
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${prefix}.${ext}
                   ${CMAKE_CURRENT_BINARY_DIR}/${prefix}.${ext} COPYONLY)
  endforeach()

  set(inputfile "${prefix}.in")

  add_test(NAME "${prefix}-serial-kokkos"
           COMMAND python ../../run_exodiff_test.py --executable "${nimble_exe}" --cli-flag "use_kokkos" --input-deck "${inputfile}" --num-ranks 1
          )

endif()
