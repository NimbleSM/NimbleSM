# ********************************************************************************************************************
#                                                      NimbleSM
#                                                VALIDATION TESTING
#
# Requirements:
# 1/ NimbleSM executables shall be named NimbleSM_XXXX with XXXX being a specific version (Serial, Kokkos, ...)
# 2/ Convention over configuration: Each test shall be stored in a separate folder. The name of the folder and the
#    files it contains must be identical.
#    Example: test/example1/example1.in example1.g example1.gold.e example1.exodiff
# 3/ All tests are run using the default python interpreter in the current environment
# 4/ A typical test contains 4 files:
#      .g       the model
#      .in      details of the simulation
#      .exodiff tolerances between gold.e file and generated .e file
#      .gold.e

# Remove below?
#if(NIMBLE_HAVE_CHARM)
#  set(charmrun ${PROJECT_BINARY_DIR}/src/charmrun)
#else()
#  set(charmrun "not available")
#endif()
#

message(STATUS "### TESTS ###")

# Defines all versions of NimbleSM to be tested
###### MISSING #if(NIMBLE_HAVE_XXXXX) to define which tests should be run
list(APPEND lstVersions
        Serial
        ArborX
#        MPI
#        Kokkos
#        Darma
#        Qthreads
#        Tpetra
        )

# Generate paths to executables
foreach(version IN LISTS lstVersions)
    set(${version}_Path ${PROJECT_BINARY_DIR}/src/NimbleSM_${version})
    message(STATUS "${version} version's path: ${${version}_Path}")
endforeach()

### Defines tests to be run for each version of NimbleSM ###

### SERIAL ###

list(APPEND lstTestsCLI_Serial
        "../run_exodiff_test.py --executable _EXECUTABLE_ --input-deck _INPUT_DECK_"
        )

list(APPEND lstTests_Serial
        brick_with_fiber_periodic
        brick_with_fibers
#        contact_entity_creation
#        cylinder_with_hole
#        plate_with_hole_rve
        quasistatic_tension
        rigid_body_motion
        simple_deformation_modes
#        single_elem_complex_displacement
#        tiny_fe_squared
#        uq_notch
#        uq_wave_in_bar
        wave_in_bar
        )

### ARBORX ###
# Currently only supports Serial version of ArborX

list(APPEND lstTestsCLI_ArborX
        "../run_exodiff_test.py --executable _EXECUTABLE_ --input-deck _INPUT_DECK_ --num-ranks 1"
        )

list(APPEND lstTests_ArborX
#        brick_with_fiber_periodic
        brick_with_fibers
        contact_entity_creation
#        cylinder_with_hole
#        notched_plate_native_hypoelastic
#        notched_plate_native_neohookean
#        plate_with_hole_rve
#        quasistatic_tension
        rigid_body_motion
        simple_deformation_modes
#        single_elem_complex_displacement
#        tiny_fe_squared
#        uq_notch
#        uq_wave_in_bar
        wave_in_bar
        )

# Build list of available tests (appending tests selected for each version of NimbleSM and removing duplicates)
foreach (version IN LISTS lstVersions)
    list(APPEND lstTests
            ${lstTests_${version}})
endforeach ()

list(REMOVE_DUPLICATES lstTests)

# NOTE failing
#add_subdirectory(cylinder_with_hole_multiblock_rve)
#add_subdirectory(plate_with_hole_rve)
#add_subdirectory(tiny_fe_squared)
#
#list(APPEND lstTests
#        brick_with_fiber_periodic
#        brick_with_fibers
#        contact_entity_creation
#        cylinder_with_hole
#        notched_plate_native_hypoelastic
#        notched_plate_native_neohookean
#        plate_with_hole_rve
#        quasistatic_tension
#        rigid_body_motion
#        simple_deformation_modes
#        single_elem_complex_displacement
#        tiny_fe_squared
#        uq_notch
#        uq_wave_in_bar
#        wave_in_bar
#        )


message(STATUS "### Copying all files necessary for testing to binary dir (extensions .in, .g, .gold.e, .exodiff)")

# Copy run_exodiff_test.py
configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/run_exodiff_test.py
        ${CMAKE_CURRENT_BINARY_DIR}/run_exodiff_test.py
        COPYONLY)

# copy all other files
foreach(test IN LISTS lstTests)
    message(STATUS " * ${test}")
    foreach(ext in g gold.e exodiff)
        configure_file(
                ${CMAKE_CURRENT_SOURCE_DIR}/${test}/${test}.${ext}
                ${CMAKE_CURRENT_BINARY_DIR}/${test}/${test}.${ext}
                COPYONLY)
    endforeach()
endforeach()


### Main loop ###
message(STATUS "## Adding tests to CMAKE: ##")
# Go through each version of NimbleSM
foreach(version IN LISTS lstVersions)
    message(STATUS " * Version: ${version}")

    # For each version, go through each test to run
    foreach(test IN LISTS lstTests_${version})

        # Determines the number of Command Line variants to run for that version of NimbleSM
        list(LENGTH lstTestsCLI_${version} count)
        math(EXPR count "${count}-1")

        # For each test, go through each Command Line variant
        foreach(i RANGE ${count})
            # Set the Command Line, test name and working directory for the current test
            list(GET lstTestsCLI_${version} ${i} cli)
            set(test_name ${test}_${version}_${i})
            set(working_dir ${CMAKE_CURRENT_BINARY_DIR}/${test})

            # Expand the command line by replacing "_EXECUTABLE_" and "_INPUT_DECK_" with variables
            string(REPLACE _EXECUTABLE_
            ${${version}_Path} cli
            ${cli}
            )

            string(REPLACE _INPUT_DECK_
            ${test}.in cli
            ${cli}
            )

            message(STATUS "  - ${test_name}: python ${cli} (work.dir: ${working_dir})")

            # Change cli from string to list (list required by add_test)
            separate_arguments(cli)

            # Add the test
            add_test(NAME ${test_name}
                    COMMAND python ${cli}
                    WORKING_DIRECTORY ${working_dir}
                    )

        endforeach() # CLI variant
    endforeach() # Test
endforeach() # NimbleSM Version
