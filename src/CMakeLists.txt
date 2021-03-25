message("src/CMakeLists.txt")

target_sources(nimble PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}/nimble_data_manager.cc
  ${CMAKE_CURRENT_LIST_DIR}/nimble_data_utils.cc
  ${CMAKE_CURRENT_LIST_DIR}/nimble_mesh_utils.cc
  ${CMAKE_CURRENT_LIST_DIR}/nimble_timing_utils.cc
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
  ${CMAKE_CURRENT_LIST_DIR}/nimble_material_factory_base.cc
  ${CMAKE_CURRENT_LIST_DIR}/nimble_model_data.cc
  ${CMAKE_CURRENT_LIST_DIR}/nimble_model_data_base.cc
  ${CMAKE_CURRENT_LIST_DIR}/nimble_expression_parser.cc
  ${CMAKE_CURRENT_LIST_DIR}/nimble_linear_solver.cc
  ${CMAKE_CURRENT_LIST_DIR}/nimble_contact_entity.cc
  ${CMAKE_CURRENT_LIST_DIR}/nimble_contact_manager.cc)

set(NIMBLE_PUBLIC_HEADERS
  ${CMAKE_CURRENT_LIST_DIR}/nimble_model_data_base.h
  ${CMAKE_CURRENT_LIST_DIR}/nimble_data_manager.h
  ${CMAKE_CURRENT_LIST_DIR}/nimble_data_utils.h
  ${CMAKE_CURRENT_LIST_DIR}/nimble_mesh_utils.h
  ${CMAKE_CURRENT_LIST_DIR}/nimble_timing_utils.h
  ${CMAKE_CURRENT_LIST_DIR}/nimble_profiling_timer.h
  ${CMAKE_CURRENT_LIST_DIR}/nimble_tpetra_utils.h
  ${CMAKE_CURRENT_LIST_DIR}/nimble_block.h
  ${CMAKE_CURRENT_LIST_DIR}/nimble_block_base.h
  ${CMAKE_CURRENT_LIST_DIR}/nimble_element.h
  ${CMAKE_CURRENT_LIST_DIR}/nimble_exodus_output.h
  ${CMAKE_CURRENT_LIST_DIR}/nimble_boundary_condition.h
  ${CMAKE_CURRENT_LIST_DIR}/nimble_parser.h
  ${CMAKE_CURRENT_LIST_DIR}/nimble_parser_util.h
  ${CMAKE_CURRENT_LIST_DIR}/nimble_boundary_condition_manager.h
  ${CMAKE_CURRENT_LIST_DIR}/nimble_genesis_mesh.h
  ${CMAKE_CURRENT_LIST_DIR}/nimble_material.h
  ${CMAKE_CURRENT_LIST_DIR}/nimble_material_factory.h
  ${CMAKE_CURRENT_LIST_DIR}/nimble_material_factory_base.h
  ${CMAKE_CURRENT_LIST_DIR}/nimble_model_data.h
  ${CMAKE_CURRENT_LIST_DIR}/nimble_block_material_interface_base.h
  ${CMAKE_CURRENT_LIST_DIR}/nimble_block_material_interface_factory_base.h
  ${CMAKE_CURRENT_LIST_DIR}/nimble_expression_parser.h
  ${CMAKE_CURRENT_LIST_DIR}/nimble_linear_solver.h
  ${CMAKE_CURRENT_LIST_DIR}/nimble_contact_interface.h
  ${CMAKE_CURRENT_LIST_DIR}/nimble_contact_entity.h
  ${CMAKE_CURRENT_LIST_DIR}/nimble_contact_manager.h
  ${CMAKE_CURRENT_LIST_DIR}/nimble_utils.h
  ${CMAKE_CURRENT_LIST_DIR}/nimble_vector_communicator.h
  ${CMAKE_CURRENT_LIST_DIR}/nimble_version.h
  ${CMAKE_CURRENT_LIST_DIR}/nimble_view.h
  ${CMAKE_CURRENT_LIST_DIR}/nimble_mpi.h)

if(NIMBLE_HAVE_KOKKOS)
  set(NIMBLE_PUBLIC_HEADERS ${NIMBLE_PUBLIC_HEADERS}
    ${CMAKE_CURRENT_LIST_DIR}/nimble_kokkos_defs.h
    ${CMAKE_CURRENT_LIST_DIR}/nimble_kokkos_block_material_interface.h
    ${CMAKE_CURRENT_LIST_DIR}/nimble_kokkos_block_material_interface_factory.h
    ${CMAKE_CURRENT_LIST_DIR}/nimble_kokkos_contact_defs.h)
endif()

if(NIMBLE_HAVE_MPI)
  set(NIMBLE_PUBLIC_HEADERS ${NIMBLE_PUBLIC_HEADERS}
    ${CMAKE_CURRENT_LIST_DIR}/nimble.mpi.mpicontext.h
    ${CMAKE_CURRENT_LIST_DIR}/nimble.mpi.rank_clique_reducer.h
    ${CMAKE_CURRENT_LIST_DIR}/nimble.mpi.reduction.h
    ${CMAKE_CURRENT_LIST_DIR}/nimble.mpi.reduction_utils.h
    ${CMAKE_CURRENT_LIST_DIR}/nimble.mpi.serialization.h
    ${CMAKE_CURRENT_LIST_DIR}/nimble.quanta.h
    ${CMAKE_CURRENT_LIST_DIR}/nimble.quanta.arrayview.h
    ${CMAKE_CURRENT_LIST_DIR}/nimble.quanta.stopwatch.h)
endif()

add_subdirectory(contact)

IF(NOT NIMBLE_HAVE_KOKKOS)
  target_sources(nimble PRIVATE ${CMAKE_CURRENT_LIST_DIR}/nimble_rve.cc)
  set(NIMBLE_PUBLIC_HEADERS ${NIMBLE_PUBLIC_HEADERS} ${CMAKE_CURRENT_LIST_DIR}/nimble_rve.h)
ENDIF()

IF(NIMBLE_HAVE_UQ)
  set(NIMBLE_PUBLIC_HEADERS ${NIMBLE_PUBLIC_HEADERS}
          ${CMAKE_CURRENT_LIST_DIR}/nimble_uq_model_data.h
          ${CMAKE_CURRENT_LIST_DIR}/nimble_uq.h )
  #
  target_sources(nimble PRIVATE
   ${CMAKE_CURRENT_LIST_DIR}/nimble_uq_model_data.cc
   ${CMAKE_CURRENT_LIST_DIR}/nimble_uq.cc
  )
ENDIF()

if (NIMBLE_HAVE_TRILINOS)
#  add_executable(NimbleSM_Tpetra nimble_tpetra.cc)
  add_executable(NimbleSM_Tpetra nimble_mpi_main.cc)
  target_link_libraries(NimbleSM_Tpetra PRIVATE nimble::nimble)
endif()

# UH -- Need to rename the executable without the MPI keyword
target_sources(nimble PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}/nimble.mpi.reduction_utils.cc
  ${CMAKE_CURRENT_LIST_DIR}/nimble.mpi.reduction.cc
  ${CMAKE_CURRENT_LIST_DIR}/nimble_mpi.cc
  )
add_executable(NimbleSM_MPI nimble_mpi_main.cc)
target_link_libraries(NimbleSM_MPI PRIVATE nimble::nimble)

IF(NIMBLE_HAVE_KOKKOS)
  target_sources(nimble PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}/nimble_kokkos_material_factory.cc
    ${CMAKE_CURRENT_LIST_DIR}/nimble_kokkos_model_data.cc
    ${CMAKE_CURRENT_LIST_DIR}/nimble_exodus_output_manager.cc
    ${CMAKE_CURRENT_LIST_DIR}/nimble_kokkos_block.cc
    ${CMAKE_CURRENT_LIST_DIR}/nimble_kokkos_block_material_interface.cc
    ${CMAKE_CURRENT_LIST_DIR}/nimble_kokkos_block_material_interface_factory.cc
  )

  set(NIMBLE_PUBLIC_HEADERS ${NIMBLE_PUBLIC_HEADERS}
    ${CMAKE_CURRENT_LIST_DIR}/nimble_kokkos_material_factory.h
    ${CMAKE_CURRENT_LIST_DIR}/nimble_kokkos_model_data.h
    ${CMAKE_CURRENT_LIST_DIR}/nimble_exodus_output_manager.h
    ${CMAKE_CURRENT_LIST_DIR}/nimble_kokkos_block.h
  )

  IF (NIMBLE_HAVE_ARBORX)
    message("Executable = NimbleSM_ArborX")
    add_executable(NimbleSM_ArborX nimble_mpi_main.cc)
    target_link_libraries(NimbleSM_ArborX PRIVATE nimble::nimble ArborX::ArborX)
  ELSE()
    message("Executable = NimbleSM_Kokkos")
    add_executable(NimbleSM_Kokkos nimble_mpi_main.cc)
    if (NIMBLE_HAVE_TRILINOS)
      target_link_libraries(NimbleSM_Kokkos nimble::nimble MPI::MPI_CXX)
    ELSE()
      if (NIMBLE_HAVE_MPI)
        target_link_libraries(NimbleSM_Kokkos nimble::nimble MPI::MPI_CXX)
      else()
        target_link_libraries(NimbleSM_Kokkos nimble::nimble)
      endif()
    ENDIF()
  ENDIF()
ENDIF()

set_target_properties(nimble PROPERTIES PUBLIC_HEADER "${NIMBLE_PUBLIC_HEADERS}")
