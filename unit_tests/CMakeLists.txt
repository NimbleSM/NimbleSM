include(CTest)

find_package(GTest REQUIRED)
include(GoogleTest)

set(NIMBLE_UNIT_SOURCES test_nimble_material_params.cc nimble_unit_main.cc)

if (NIMBLE_HAVE_KOKKOS)
  set(NIMBLE_UNIT_SOURCES ${NIMBLE_UNIT_SOURCES} test_nimble_kokkos_material.cc)
endif()

add_executable(NimbleSM_Unit ${NIMBLE_UNIT_SOURCES})

gtest_add_tests(NimbleSM_Unit SOURCES ${NIMBLE_UNIT_SOURCES})

# target_include_directories(NimbleSM_Unit PRIVATE ${GTEST_INCLUDE_DIR})
target_link_libraries(NimbleSM_Unit nimble::nimble GTest::GTest)
if (NIMBLE_HAVE_TRILINOS)
  target_link_libraries(NimbleSM_Unit Trilinos::Trilinos MPI::MPI_CXX)
else()
  if (NIMBLE_HAVE_MPI)
    target_link_libraries(NimbleSM_Unit Kokkos::kokkos MPI::MPI_CXX)
  else()
    target_link_libraries(NimbleSM_Unit Kokkos::kokkos)
  endif()
endif()
