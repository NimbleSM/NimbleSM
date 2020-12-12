
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/uq_notch.in
               ${CMAKE_CURRENT_BINARY_DIR}/uq_notch.in COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/uq_notch.g
               ${CMAKE_CURRENT_BINARY_DIR}/uq_notch.g COPYONLY)
#configure_file(${CMAKE_CURRENT_SOURCE_DIR}/uq_notch.g.3.0
#               ${CMAKE_CURRENT_BINARY_DIR}/uq_notch.g.3.0 COPYONLY)
#configure_file(${CMAKE_CURRENT_SOURCE_DIR}/uq_notch.g.3.1
#               ${CMAKE_CURRENT_BINARY_DIR}/uq_notch.g.3.1 COPYONLY)
#configure_file(${CMAKE_CURRENT_SOURCE_DIR}/uq_notch.g.3.2
#               ${CMAKE_CURRENT_BINARY_DIR}/uq_notch.g.3.2 COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/uq_notch.gold.e
               ${CMAKE_CURRENT_BINARY_DIR}/uq_notch.gold.e COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/uq_notch.exodiff
               ${CMAKE_CURRENT_BINARY_DIR}/uq_notch.exodiff COPYONLY)

add_test(uq_notch_serial python ../run_exodiff_test.py --executable "${NimbleSM_Serial}" --input-deck uq_notch.in)
