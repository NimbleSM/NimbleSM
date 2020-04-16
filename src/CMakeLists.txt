
target_sources(nimble PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}/nimble_data_manager.cc
    ${CMAKE_CURRENT_LIST_DIR}/nimble_data_utils.cc
    ${CMAKE_CURRENT_LIST_DIR}/nimble_mesh_utils.cc
    ${CMAKE_CURRENT_LIST_DIR}/nimble_tpetra_utils.cc
    ${CMAKE_CURRENT_LIST_DIR}/nimble_block.cc
    ${CMAKE_CURRENT_LIST_DIR}/nimble_element.cc
    ${CMAKE_CURRENT_LIST_DIR}/nimble_exodus_output.cc
    ${CMAKE_CURRENT_LIST_DIR}/nimble_boundary_condition.cc
    ${CMAKE_CURRENT_LIST_DIR}/nimble_parser.cc
    ${CMAKE_CURRENT_LIST_DIR}/nimble_boundary_condition_manager.cc
    ${CMAKE_CURRENT_LIST_DIR}/nimble_genesis_mesh.cc
    ${CMAKE_CURRENT_LIST_DIR}/nimble_material.cc
    ${CMAKE_CURRENT_LIST_DIR}/nimble_material_factory.cc
    ${CMAKE_CURRENT_LIST_DIR}/nimble_material_factory_util.cc
    ${CMAKE_CURRENT_LIST_DIR}/nimble_expression_parser.cc
    ${CMAKE_CURRENT_LIST_DIR}/nimble_linear_solver.cc
    ${CMAKE_CURRENT_LIST_DIR}/nimble_contact_entity.cc
    ${CMAKE_CURRENT_LIST_DIR}/nimble_contact_manager.cc
    )

IF(NOT NIMBLE_HAVE_KOKKOS)
  target_sources(nimble PRIVATE ${CMAKE_CURRENT_LIST_DIR}/nimble_rve.cc)
ENDIF()

IF(NIMBLE_HAVE_UQ)
  target_sources(nimble PRIVATE ${CMAKE_CURRENT_LIST_DIR}/nimble_uq.cc)
ENDIF()

add_executable(NimbleSM_Serial nimble_serial.cc)
target_link_libraries(NimbleSM_Serial PRIVATE nimble::nimble)
if (NIMBLE_HAVE_EXTRAS)
  target_link_libraries(NimbleSM_Serial PRIVATE nimble::extras)
endif()


if (NIMBLE_HAVE_TRILINOS)
  add_executable(NimbleSM_Tpetra nimble_tpetra.cc)
  target_link_libraries(NimbleSM_Tpetra PRIVATE nimble::nimble)
  if (NIMBLE_HAVE_EXTRAS)
    target_link_libraries(NimbleSM_Tpetra PRIVATE nimble::extras)
  endif()
endif()

IF(NIMBLE_HAVE_MPI)
  add_executable(NimbleSM_MPI
      nimble_mpi.cc
      nimble.mpi.reduction_utils.cc
      nimble.mpi.reduction.cc
      )
  target_link_libraries(NimbleSM_MPI PRIVATE nimble::nimble )
  if (NIMBLE_HAVE_EXTRAS)
    target_link_libraries(NimbleSM_MPI PRIVATE nimble::extras)
  endif()
ENDIF()

IF(NIMBLE_HAVE_QTHREADS)
  add_executable(NimbleSM_Qthreads nimble_qthread.cc ${NIMBLE_SOURCES})
  target_link_libraries(NimbleSM_Qthreads ${LAME_LIBS} ${Exodus_LIBS} ${bvh_LIBS} ${Qthreads_LIBS})
ENDIF()

IF(NIMBLE_HAVE_KOKKOS)
  target_sources(nimble PRIVATE ${CMAKE_CURRENT_LIST_DIR}/nimble_kokkos_material_factory.cc
                 ${CMAKE_CURRENT_LIST_DIR}/nimble_kokkos_data_manager.cc
                 ${CMAKE_CURRENT_LIST_DIR}/nimble_exodus_output_manager.cc
                 ${CMAKE_CURRENT_LIST_DIR}/nimble_kokkos_block.cc
)
  add_executable(NimbleSM_Kokkos nimble_kokkos.cc nimble.mpi.reduction.cc nimble.mpi.reduction_utils.cc)
  target_link_libraries(NimbleSM_Kokkos nimble::nimble Trilinos::Trilinos MPI::MPI_CXX)
  if (NIMBLE_HAVE_EXTRAS)
    target_link_libraries(NimbleSM_Kokkos nimble::extras)
  endif()
ENDIF()
