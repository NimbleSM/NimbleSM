
set(prefix "brick_with_fiber_periodic")

if(NIMBLE_HAVE_TRILINOS)

  foreach (ext "in" "g" "gold.e" "exodiff")
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${prefix}.${ext}
                   ${CMAKE_CURRENT_BINARY_DIR}/${prefix}.${ext} COPYONLY)
  endforeach()

  set(inputfile "${prefix}.in")

  add_test(NAME "${prefix}-serial-tpetra"
           COMMAND python ../../run_test.py --executable "${nimble_exe}" --cli-flag "use_tpetra" --input-deck "${inputfile}" --num-ranks 1
          )

endif()

