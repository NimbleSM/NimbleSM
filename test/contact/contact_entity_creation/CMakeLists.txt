
set(prefix "contact_entity_creation")

set(inputfile "${prefix}.in")

if (NIMBLE_HAVE_KOKKOS)

  foreach (ext "in" "g" "gold.e" "exodiff")
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${prefix}.${ext}
                   ${CMAKE_CURRENT_BINARY_DIR}/${prefix}.${ext} COPYONLY)
  endforeach()

  add_test(NAME "${prefix}-serial-kokkos"
           COMMAND python ../../run_exodiff_test.py --executable "${nimble_exe}" --cli-flag "use_kokkos" --input-deck "${inputfile}" --num-ranks 1
          )

  if (NIMBLE_HAVE_MPI)

    foreach (ext "g.2.0" "g.2.1" "g.4.0" "g.4.1" "g.4.2" "g.4.3")
      configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${prefix}.${ext}
                     ${CMAKE_CURRENT_BINARY_DIR}/${prefix}.${ext} COPYONLY)
    endforeach()

    foreach (nrank 2 4)

      add_test(NAME "${prefix}-np${nrank}-kokkos"
               COMMAND python ../../run_exodiff_test.py --executable "${nimble_exe}" --cli-flag "use_kokkos" --input-deck "${inputfile}" --num-ranks ${nrank}
              )
    endforeach()

  endif()

endif()

