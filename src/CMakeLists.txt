
target_sources(nimble PRIVATE
    nimble_data_manager.cc
    nimble_data_utils.cc
    nimble_mesh_utils.cc
    nimble_tpetra_utils.cc
    nimble_block.cc
    nimble_element.cc
    nimble_exodus_output.cc
    nimble_boundary_condition.cc
    nimble_parser.cc
    nimble_boundary_condition_manager.cc
    nimble_genesis_mesh.cc
    nimble_material.cc
    nimble_material_factory.cc
    nimble_material_factory_util.cc
    nimble_expression_parser.cc
    nimble_linear_solver.cc
    nimble_contact_entity.cc
    nimble_contact_manager.cc
    )

IF(NOT NIMBLE_HAVE_KOKKOS)
  target_sources(nimble PRIVATE nimble_rve.cc)
ENDIF()

set(NIMBLE_KOKKOS_SOURCES
  nimble_kokkos_data_manager.cc
  nimble_kokkos_block.cc
  nimble_kokkos_material_factory.cc
  nimble_exodus_output_manager.cc
)

add_executable(NimbleSM_Serial nimble_serial.cc)
target_link_libraries(NimbleSM_Serial PRIVATE nimble::nimble)

set(NIMBLE_MAINS
  nimble_serial.cc nimble_tpetra.cc nimble_mpi.cc nimble_qthread.cc nimble_darma.cc nimble_kokkos.cc
)

if (NIMBLE_NVIDIA_BUILD)
  set_source_files_properties(${NIMBLE_SOURCES} PROPERTIES LANGUAGE CUDA)
  set_source_files_properties(${NIMBLE_DARMA_SOURCES} PROPERTIES LANGUAGE CUDA)
  set_source_files_properties(${NIMBLE_KOKKOS_SOURCES} PROPERTIES LANGUAGE CUDA)
  set_source_files_properties(${NIMBLE_TPETRA_SOURCES} PROPERTIES LANGUAGE CUDA)
  set_source_files_properties(${NIMBLE_MPI_SOURCES} PROPERTIES LANGUAGE CUDA)
  set_source_files_properties(${NIMBLE_MAINS} PROPERTIES LANGUAGE CUDA)
  set_source_files_properties(${NIMBLE_EXTRAS_SOURCES} PROPERTIES LANGUAGE CUDA)
endif()

set(NIMBLE_SOURCES ${NIMBLE_SOURCES} ${NIMBLE_EXTRAS_SOURCES})

if (NIMBLE_HAVE_MPI OR NIMBLE_HAVE_QTHREADS OR NIMBLE_HAVE_DARMA OR NIMBLE_HAVE_KOKKOS)
  set(NIMBLE_SOURCES ${NIMBLE_SOURCES} ${NIMBLE_MPI_SOURCES})
endif()

if (NIMBLE_HAVE_DARMA)
  set(NIMBLE_SOURCES ${NIMBLE_SOURCES} ${NIMBLE_DARMA_SOURCES})
endif()

if(NIMBLE_HAVE_KOKKOS)
  set(NIMBLE_SOURCES ${NIMBLE_SOURCES} ${NIMBLE_KOKKOS_SOURCES})
endif()

add_executable(NimbleSM_Serial nimble_serial.cc ${NIMBLE_SOURCES})
target_link_libraries(NimbleSM_Serial ${NIMBLE_EXTRAS_LIBS} ${NIMBLE_KOKKOS_LIBS} ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES} ${Exodus_LIBS} ${bvh_LIBS} ${MPI_CXX_LIBRARIES})
install(TARGETS NimbleSM_Serial DESTINATION bin)

IF(NIMBLE_HAVE_TRILINOS)
  add_executable(NimbleSM_Tpetra nimble_tpetra.cc ${NIMBLE_SOURCES})
  target_link_libraries(NimbleSM_Tpetra ${NIMBLE_EXTRAS_LIBS} ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES} ${Exodus_LIBS} ${bvh_LIBS} ${MPI_CXX_LIBRARIES})
ENDIF()

IF(NIMBLE_HAVE_MPI)
  add_executable(NimbleSM_MPI nimble_mpi.cc ${NIMBLE_SOURCES})
  target_link_libraries(NimbleSM_MPI ${NIMBLE_EXTRAS_LIBS} ${NIMBLE_KOKKOS_LIBS} ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES} ${Exodus_LIBS} ${bvh_LIBS} ${MPI_CXX_LIBRARIES})
ENDIF()

IF(NIMBLE_HAVE_QTHREADS)
  add_executable(NimbleSM_Qthreads nimble_qthread.cc ${NIMBLE_SOURCES})
  target_link_libraries(NimbleSM_Qthreads ${LAME_LIBS} ${Exodus_LIBS} ${bvh_LIBS} ${Qthreads_LIBS})
ENDIF()

IF(NIMBLE_HAVE_DARMA)
  add_executable(NimbleSM_Darma nimble_darma.cc ${NIMBLE_SOURCES})
  IF(NIMBLE_HAVE_TRILINOS)
    target_link_libraries(NimbleSM_Darma ${NIMBLE_EXTRAS_LIBS} ${Exodus_LIBS} ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES} ${DARMALIBRARIES} ${bvh_LIBS})
  ELSE()
    target_link_libraries(NimbleSM_Darma ${NIMBLE_EXTRAS_LIBS} ${Exodus_LIBS} ${DARMA_LIBRARIES} ${bvh_LIBS})
    target_link_libraries(NimbleSM_Serial ${NIMBLE_EXTRAS_LIBS} ${Exodus_LIBS} ${DARMA_LIBRARIES} ${bvh_LIBS})
  ENDIF()
  install(TARGETS NimbleSM_Darma DESTINATION bin)
ENDIF()

IF(NIMBLE_HAVE_KOKKOS)
  add_executable(NimbleSM_Kokkos nimble_kokkos.cc ${NIMBLE_SOURCES})
  target_link_libraries(NimbleSM_Kokkos ${NIMBLE_EXTRAS_LIBS} ${NIMBLE_KOKKOS_LIBS} ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES} ${Exodus_LIBS} ${bvh_LIBS} ${MPI_CXX_LIBRARIES})
ENDIF()
