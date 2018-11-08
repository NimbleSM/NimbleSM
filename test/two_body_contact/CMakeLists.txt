
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/two_body_contact.in
               ${CMAKE_CURRENT_BINARY_DIR}/two_body_contact.in COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/two_body_contact.g
               ${CMAKE_CURRENT_BINARY_DIR}/two_body_contact.g COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/two_body_contact.gold.e
               ${CMAKE_CURRENT_BINARY_DIR}/two_body_contact.gold.e COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/two_body_contact.exodiff
               ${CMAKE_CURRENT_BINARY_DIR}/two_body_contact.exodiff COPYONLY)

if(NIMBLE_HAVE_EXTRAS)
  add_test(two_body_contact_serial python ../run_exodiff_test.py --executable "${NimbleSM_Serial}" --input-deck two_body_contact.in)
endif()
