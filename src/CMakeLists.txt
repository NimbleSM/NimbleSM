
set(NIMBLE_SOURCES
  nimble_data_manager.cc
  nimble_data_utils.cc
  nimble_mesh_utils.cc
  nimble_tpetra_utils.cc
  nimble.mpi.utils.cc
  nimble_block.cc
  nimble_element.cc
  nimble_exodus_output.cc
  nimble_boundary_condition.cc
  nimble_parser.cc
  nimble_boundary_condition_manager.cc
  nimble_genesis_mesh.cc
  nimble_material.cc
  nimble_expression_parser.cc
  nimble_linear_solver.cc
  nimble_contact.cc
)

IF(NOT NIMBLE_HAVE_KOKKOS)
  set(NIMBLE_SOURCES ${NIMBLE_SOURCES} nimble_rve.cc)
ENDIF()

set(NIMBLE_KOKKOS_SOURCES
  nimble_kokkos_data_manager.cc
  nimble_kokkos_block.cc
  nimble_exodus_output_manager.cc
)

set(NIMBLE_DARMA_SOURCES
  nimble_darma_utils.cc
)

add_executable(NimbleSM_Serial nimble_serial.cc ${NIMBLE_SOURCES} ${NIMBLE_EXTRAS_SOURCES})
target_link_libraries(NimbleSM_Serial ${NIMBLE_KOKKOS_LIBS} ${NIMBLE_EXTRAS_LIBS} ${Trilinos_TPL_LIBRARIES} ${Trilinos_LIBRARIES} ${Exodus_LIBS} ${bvh_LIBS})
install(TARGETS NimbleSM_Serial DESTINATION bin)

IF(NIMBLE_HAVE_TRILINOS)
  add_executable(NimbleSM_Tpetra nimble_tpetra.cc ${NIMBLE_SOURCES} ${NIMBLE_EXTRAS_SOURCES})
  target_link_libraries(NimbleSM_Tpetra ${NIMBLE_EXTRAS_LIBS} ${Trilinos_TPL_LIBRARIES} ${Trilinos_LIBRARIES} ${bvh_LIBS})
ENDIF()

IF(NIMBLE_HAVE_MPI)
  add_executable(NimbleSM_MPI nimble_mpi.cc ${NIMBLE_SOURCES} ${NIMBLE_EXTRAS_SOURCES})
  target_link_libraries(NimbleSM_MPI ${NIMBLE_KOKKOS_LIBS} ${NIMBLE_EXTRAS_LIBS} ${Exodus_LIBS} ${bvh_LIBS})
ENDIF()

IF(NIMBLE_HAVE_QTHREADS)
  add_executable(NimbleSM_Qthreads nimble_qthread.cc ${NIMBLE_SOURCES} ${NIMBLE_EXTRAS_SOURCES})
  target_link_libraries(NimbleSM_Qthreads ${Qthreads_LIBS} ${LAME_LIBS} ${Exodus_LIBS} ${bvh_LIBS})
ENDIF()

IF(NIMBLE_HAVE_DARMA)
  add_executable(NimbleSM_Darma nimble_darma.cc ${NIMBLE_SOURCES} ${NIMBLE_EXTRAS_SOURCES} ${NIMBLE_DARMA_SOURCES})
  IF(NIMBLE_HAVE_TRILINOS)
    target_link_libraries(NimbleSM_Darma ${NIMBLE_EXTRAS_LIBS} ${Trilinos_TPL_LIBRARIES} ${Trilinos_LIBRARIES} ${DARMALIBRARIES} ${bvh_LIBS})
  ELSE()
    target_link_libraries(NimbleSM_Darma ${NIMBLE_EXTRAS_LIBS} ${Exodus_LIBS} ${DARMA_LIBRARIES} ${bvh_LIBS})
    target_link_libraries(NimbleSM_Serial ${NIMBLE_EXTRAS_LIBS} ${Exodus_LIBS} ${DARMA_LIBRARIES} ${bvh_LIBS})
  ENDIF()
  install(TARGETS NimbleSM_Darma DESTINATION bin)
ENDIF()

IF(NIMBLE_HAVE_KOKKOS)
  add_executable(KokkosSandbox nimble_kokkos_sandbox.cc ${NIMBLE_SOURCES} ${NIMBLE_EXTRAS_SOURCES} ${NIMBLE_KOKKOS_SOURCES} )
  target_link_libraries(KokkosSandbox ${NIMBLE_EXTRAS_LIBS} ${NIMBLE_KOKKOS_LIBS} ${Exodus_LIBS})
ENDIF()
