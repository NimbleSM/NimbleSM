message("### UNIT TESTS CONFIGURATION STARTS ###")

message("- Necessary cmake includes (CTest, GoogleTest, FetchContent)")
include(CTest)
include(GoogleTest)
include(FetchContent)

set(NIMBLE_USE_LOCAL_GTEST OFF CACHE BOOL "Whether to use local GTest or fetch from Github")

if (NIMBLE_USE_LOCAL_GTEST)
  message("- Using local Google Test")
  find_package(GTest REQUIRED)
else()
  # Fetch Google Test
  set(GOOGLETEST_VERSION 1.10.0)
  set(GOOGLETEST_GIT_URL https://github.com/google/googletest.git)

  message("- Using FetchContent to get Google Test version ${GOOGLETEST_VERSION} from ${GOOGLETEST_GIT_URL}")

  FetchContent_Declare(
          googletest
          GIT_REPOSITORY ${GOOGLETEST_GIT_URL}
          GIT_TAG        release-${GOOGLETEST_VERSION}
  )

  FetchContent_MakeAvailable(googletest)
endif()

message(" * gtest sources dir: ${gtest_SOURCE_DIR}")
message(" * gtest binary dir:  ${gtest_BINARY_DIR}")

# Configure tests
message("- Configure tests")

set(NIMBLE_UNIT_SOURCES
        nimble_unit_main.cc
        projection_node_to_face.cc
        test_nimble_material_params.cc
        test_nimble_utils.cc
        )

if (NIMBLE_HAVE_KOKKOS)
    set(NIMBLE_UNIT_SOURCES ${NIMBLE_UNIT_SOURCES} test_nimble_kokkos_material.cc)
endif()

message(" * The following test files will be added: ${NIMBLE_UNIT_SOURCES}")
message(" * Call to gtest_add_tests")
gtest_add_tests(NimbleSM_Unit SOURCES ${NIMBLE_UNIT_SOURCES})


# Configure Executable
message("- Configure new executable NimbleSM_Unit")
message(" * Add Executable using the following source files:  ${NIMBLE_UNIT_SOURCES}")
add_executable(NimbleSM_Unit ${NIMBLE_UNIT_SOURCES})
message(" * Link nimble library")
target_link_libraries(NimbleSM_Unit nimble::nimble)
message(" * Link gtest libraries")
if (NIMBLE_USE_LOCAL_GTEST)
  target_include_directories(NimbleSM_Unit PRIVATE ${GTEST_INCLUDE_DIR})
  target_link_libraries(NimbleSM_Unit ${GTEST_LIBRARY})
else()
  target_link_libraries(NimbleSM_Unit gtest gtest_main)
endif()

message(" * Link other libraries")

if (NIMBLE_HAVE_TRILINOS)
  target_link_libraries(NimbleSM_Unit Trilinos::Trilinos)
else()

  if (NIMBLE_HAVE_MPI)
    target_link_libraries(NimbleSM_Unit MPI::MPI_CXX)
  endif()

  if(NIMBLE_HAVE_KOKKOS)
    target_link_libraries(NimbleSM_Unit Kokkos::kokkos)
  endif()
endif()

message("### UNIT TESTS CONFIGURATION COMPLETED ###")
