
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/contact_entity_creation.in
               ${CMAKE_CURRENT_BINARY_DIR}/contact_entity_creation.in COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/contact_entity_creation.g
               ${CMAKE_CURRENT_BINARY_DIR}/contact_entity_creation.g COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/contact_entity_creation.g.2.0
               ${CMAKE_CURRENT_BINARY_DIR}/contact_entity_creation.g.2.0 COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/contact_entity_creation.g.2.1
               ${CMAKE_CURRENT_BINARY_DIR}/contact_entity_creation.g.2.1 COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/contact_entity_creation.g.3.0
               ${CMAKE_CURRENT_BINARY_DIR}/contact_entity_creation.g.3.0 COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/contact_entity_creation.g.3.1
               ${CMAKE_CURRENT_BINARY_DIR}/contact_entity_creation.g.3.1 COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/contact_entity_creation.g.3.2
               ${CMAKE_CURRENT_BINARY_DIR}/contact_entity_creation.g.3.2 COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/contact_entity_creation.g.4.0
               ${CMAKE_CURRENT_BINARY_DIR}/contact_entity_creation.g.4.0 COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/contact_entity_creation.g.4.1
               ${CMAKE_CURRENT_BINARY_DIR}/contact_entity_creation.g.4.1 COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/contact_entity_creation.g.4.2
               ${CMAKE_CURRENT_BINARY_DIR}/contact_entity_creation.g.4.2 COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/contact_entity_creation.g.4.3
               ${CMAKE_CURRENT_BINARY_DIR}/contact_entity_creation.g.4.3 COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/contact_entity_creation.gold.e
               ${CMAKE_CURRENT_BINARY_DIR}/contact_entity_creation.gold.e COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/contact_entity_creation.exodiff
               ${CMAKE_CURRENT_BINARY_DIR}/contact_entity_creation.exodiff COPYONLY)

if(NIMBLE_HAVE_KOKKOS)
  add_test(contact_entity_creation_kokkos_np1 python ../run_exodiff_test.py --executable "${NimbleSM_Kokkos}" --input-deck contact_entity_creation.in --num-ranks 1)
  add_test(contact_entity_creation_kokkos_np2 python ../run_exodiff_test.py --executable "${NimbleSM_Kokkos}" --input-deck contact_entity_creation.in --num-ranks 2)
  add_test(contact_entity_creation_kokkos_np3 python ../run_exodiff_test.py --executable "${NimbleSM_Kokkos}" --input-deck contact_entity_creation.in --num-ranks 3)
  add_test(contact_entity_creation_kokkos_np4 python ../run_exodiff_test.py --executable "${NimbleSM_Kokkos}" --input-deck contact_entity_creation.in --num-ranks 4)
endif()
