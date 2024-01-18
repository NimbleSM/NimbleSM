
set(prefix "cubes_contact")

foreach (ext "in" "g" "gold.e" "exodiff")
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${prefix}.${ext}
                 ${CMAKE_CURRENT_BINARY_DIR}/${prefix}.${ext} COPYONLY)
endforeach()

set(inputfile "${prefix}.in")

if (NIMBLE_HAVE_ARBORX)
  #
  # When ArborX is set, Kokkos is also present
  #

  add_test(NAME "${prefix}-serial-kokkos"
           COMMAND Python::Interpreter ../../run_exodiff_test.py --executable "${nimble_exe}" --input-deck "${inputfile}" --num-ranks 1
          )

endif()

if (NIMBLE_HAVE_BVH)
  #
  # When BVH is set, VT is also present
  #
  add_test(NAME "${prefix}-serial-bvh"
           COMMAND Python::Interpreter ../../run_exodiff_test.py --executable "${nimble_exe}" --cli-flag "use_vt" --input-deck "${inputfile}" --num-ranks 1
          )
endif()

