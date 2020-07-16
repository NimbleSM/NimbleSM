include(CTest)

find_package(GTest REQUIRED)
include(GoogleTest)

set(NIMBLE_UNIT_SOURCES test_nimble_material_params.cc nimble_unit_main.cc)

add_executable(NimbleSM_Unit ${NIMBLE_UNIT_SOURCES})

gtest_discover_tests(NimbleSM_Unit)

target_include_directories(NimbleSM_Unit PRIVATE ${GTEST_INCLUDE_DIR})
target_link_libraries(NimbleSM_Unit PRIVATE nimble::nimble ${GTEST_LIBRARY})
