
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/quasistatic_tension.in
               ${CMAKE_CURRENT_BINARY_DIR}/quasistatic_tension.in COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/quasistatic_tension.g
               ${CMAKE_CURRENT_BINARY_DIR}/quasistatic_tension.g COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/quasistatic_tension.g.4.0
               ${CMAKE_CURRENT_BINARY_DIR}/quasistatic_tension.g.4.0 COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/quasistatic_tension.g.4.1
               ${CMAKE_CURRENT_BINARY_DIR}/quasistatic_tension.g.4.1 COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/quasistatic_tension.g.4.2
               ${CMAKE_CURRENT_BINARY_DIR}/quasistatic_tension.g.4.2 COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/quasistatic_tension.g.4.3
               ${CMAKE_CURRENT_BINARY_DIR}/quasistatic_tension.g.4.3 COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/quasistatic_tension.gold.e
               ${CMAKE_CURRENT_BINARY_DIR}/quasistatic_tension.gold.e COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/quasistatic_tension.exodiff
               ${CMAKE_CURRENT_BINARY_DIR}/quasistatic_tension.exodiff COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/PostProcess.py
               ${CMAKE_CURRENT_BINARY_DIR}/PostProcess.py COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/quasistatic_tension.gnuplot
               ${CMAKE_CURRENT_BINARY_DIR}/quasistatic_tension.gnuplot COPYONLY)

if(NIMBLE_HAVE_EXTRAS)
  add_test(quasistatic_tension_J2_lame_serial python ../run_exodiff_test.py --executable "${NimbleSM_Serial}" --input-deck quasistatic_tension.in)
endif()
