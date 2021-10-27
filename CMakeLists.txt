#
# CMakeLists.txt file for NimbleSM
#

message("CMakeLists.txt - NimbleSM")

#--- Set CMake minimum version
cmake_minimum_required(VERSION 3.17)

project(NimbleSM VERSION 0.1.0 LANGUAGES CXX)

#--- Define options

option(NimbleSM_ENABLE_ARBORX "Whether to use ArborX" OFF)
option(HAVE_ARBORX "Whether to use ArborX" OFF)
if (HAVE_ARBORX)
  set(NimbleSM_ENABLE_ARBORX TRUE)
  message(DEPRECATION "Flag HAVE_ARBORX is deprecated. Use NimbleSM_ENABLE_ARBORX instead.")
endif()
#------

option(NimbleSM_ENABLE_BVH "Whether to use BVH for parallel asynchronous contact" OFF)
option(HAVE_BVH "Whether to use BVH for parallel asynchronous contact" OFF)
if (HAVE_BVH)
  set(NimbleSM_ENABLE_BVH TRUE)
  message(DEPRECATION "Flag HAVE_BVH is deprecated. Use NimbleSM_ENABLE_BVH instead.")
endif()
#------


option(NimbleSM_ENABLE_CUDA "Whether to build Nimble with CUDA" OFF)
option(NIMBLE_NVIDIA_BUILD "Whether to build Nimble with CUDA" OFF)
if (NIMBLE_NVIDIA_BUILD)
  set(NimbleSM_ENABLE_CUDA TRUE)
  message(DEPRECATION "Flag NIMBLE_NVIDIA_BUILD is deprecated. Use NimbleSM_ENABLE_CUDA instead.")
endif()
if (NimbleSM_ENABLE_CUDA)
  message(WARNING "NVIDIA build is not currently active. ")
  #add_definitions(-DNIMBLE_NVIDIA_BUILD)
endif()
#------

option(NimbleSM_ENABLE_UNIT_TESTS "Whether to build Nimble with unit testing" OFF)
option(NIMBLE_ENABLE_UNIT_TESTS "Whether to build Nimble with unit testing" OFF)
if (NIMBLE_ENABLE_UNIT_TESTS)
  set(NimbleSM_ENABLE_UNIT_TESTS TRUE)
  message(DEPRECATION "Flag NIMBLE_ENABLE_UNIT_TESTS is deprecated. Use NimbleSM_ENABLE_UNIT_TESTS instead.")
endif()
#------

option(NimbleSM_ENABLE_KOKKOS "Whether to use Kokkos" OFF)
option(HAVE_KOKKOS "Whether to use Kokkos" OFF)
if (HAVE_KOKKOS)
  set(NimbleSM_ENABLE_KOKKOS TRUE)
  message(DEPRECATION "Flag HAVE_KOKKOS is deprecated. Use NimbleSM_ENABLE_KOKKOS instead.")
endif()
#------

option(NimbleSM_ENABLE_MPI "Whether to use Pure MPI (no Trilinos)" OFF)
option(USE_PURE_MPI "Whether to use Pure MPI (no Trilinos)" OFF)
if (USE_PURE_MPI)
  set(NimbleSM_ENABLE_MPI TRUE)
  message(DEPRECATION "Flag USE_PURE_MPI is deprecated. Use NimbleSM_ENABLE_MPI instead.")
endif()
#------

option(NimbleSM_ENABLE_TRILINOS "Whether to use Trilinos" OFF)
option(HAVE_TRILINOS "Whether to use Trilinos" OFF)
if (HAVE_TRILINOS)
  set(NimbleSM_ENABLE_TRILINOS TRUE)
  message(DEPRECATION "Flag HAVE_TRILINOS is deprecated. Use NimbleSM_ENABLE_TRILINOS instead.")
endif()
#------

option(NimbleSM_ENABLE_UQ "Whether to use enamble UQ sampling" OFF)
option(HAVE_UQ "Whether to use enamble UQ sampling" OFF)
if (HAVE_UQ)
  set(NimbleSM_ENABLE_UQ TRUE)
  message(DEPRECATION "Flag HAVE_UQ is deprecated. Use NimbleSM_ENABLE_UQ instead.")
endif()
#------

option(NimbleSM_TIME_Contact "Whether to time contact" OFF)
option(TIME_CONTACT "Whether to time contact" OFF)
if (TIME_CONTACT)
  set(NimbleSM_TIME_Contact TRUE)
  message(DEPRECATION "Flag TIME_CONTACT is deprecated. Use NimbleSM_TIME_Contact instead.")
endif()
#------

include(GNUInstallDirs)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

enable_testing()

set(SRC_DIR ${CMAKE_SOURCE_DIR}/src)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

add_library(nimble)
add_library(nimble::nimble ALIAS nimble)

# Enhanced debugging option
OPTION(NimbleSM_ENABLE_DEBUG "Flag to turn on enhanced debugging info" OFF)
IF (NimbleSM_ENABLE_DEBUG)
  MESSAGE("-- Enhanced debugging        is Enabled, compiling with -DNIMBLE_DEBUG")
  SET(NIMBLE_DEBUG TRUE)
ELSE()
  SET(NIMBLE_DEBUG FALSE)
ENDIF()

target_include_directories(nimble PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

# Optional functionality for UQ
if (NimbleSM_ENABLE_UQ)
  set(NIMBLE_HAVE_UQ TRUE)
  message(STATUS "Compiling with UQ")
  target_compile_definitions(nimble PUBLIC NIMBLE_HAVE_UQ)
endif()

# Optional functionality for contact
if (NimbleSM_TIME_Contact)
  set(NIMBLE_TIME_CONTACT TRUE)
  message(STATUS "Compiling with contact timing")
  target_compile_definitions(nimble PUBLIC NIMBLE_TIME_CONTACT)
endif()

# Optional dependency for pure MPI build (no Trilinos)
if (NimbleSM_ENABLE_MPI)
  set(NIMBLE_HAVE_MPI TRUE)
  find_package(MPI REQUIRED)
  target_link_libraries(nimble PUBLIC MPI::MPI_CXX)
  target_compile_definitions(nimble PUBLIC NIMBLE_HAVE_MPI)
  message(STATUS "Compiling with MPI")
endif()

# Optional BVH dependency
if (NimbleSM_ENABLE_BVH)
  ADD_DEFINITIONS(-DNIMBLE_HAVE_BVH)
  target_compile_definitions(nimble PUBLIC -DNIMBLE_HAVE_VT)
  if (BVH_DEBUG_LEVEL)
    ADD_DEFINITIONS(-DBVH_DEBUG_LEVEL=${BVH_DEBUG_LEVEL})
  endif()
  set(NIMBLE_HAVE_BVH TRUE)
  set(bvh_DIR CACHE PATH "Path to BVH install directory")
  find_package(bvh REQUIRED)
  message(STATUS "Compiling with BVH:")
  message("       BVH_DIR: ${bvh_DIR}")
  if (NIMBLE_HAVE_MPI)
    message("       Using parallel collision detection (NIMBLE_HAVE_MPI=ON)")
  else()
    message("       Using only serial collision detection (NIMBLE_HAVE_MPI=OFF)")
  endif()
  target_link_libraries(nimble PUBLIC bvh::bvh)
else()
  message(STATUS "BVH is NOT enabled.")
  SET(NIMBLE_HAVE_BVH FALSE)
  set(bvh_LIBS)
endif()

# Optional Trilinos dependencies
if(NimbleSM_ENABLE_TRILINOS)
  option(NimbleSM_USE_TRILINOS_EXODUS "Whether to use Exodus from Trilinos" ON)
  option(USE_TRILINOS_EXODUS "Whether to use Exodus from Trilinos" ON)
  if (NOT USE_TRILINOS_EXODUS)
    set(NimbleSM_USE_TRILINOS_EXODUS FALSE)
    message(DEPRECATION "Flag USE_TRILINOS_EXODUS is deprecated. Use NimbleSM_USE_TRILINOS_EXODUS instead.")
  endif()
  # Trilinos brings in Trilinos, Kokkos, and (optionally) Exodus
  target_compile_definitions(nimble PUBLIC -DNIMBLE_HAVE_TRILINOS)
  set(NIMBLE_HAVE_TRILINOS TRUE)
  target_compile_definitions(nimble PUBLIC -DNIMBLE_HAVE_KOKKOS)
  set(NIMBLE_HAVE_KOKKOS TRUE)
  if(NimbleSM_USE_TRILINOS_EXODUS)
    target_compile_definitions(nimble PUBLIC -DNIMBLE_HAVE_EXODUS)
    set(NIMBLE_HAVE_EXODUS TRUE)
  endif()
  message(STATUS "Trilinos is enabled.")
  set(Trilinos_DIR CACHE PATH "Path to Trilinos install directory")
  set(NIMBLE_TRILINOS_COMPONENTS KokkosCore TpetraCore)
  if(NimbleSM_USE_TRILINOS_EXODUS)
    set(NIMBLE_TRILINOS_COMPONENTS ${NIMBLE_TRILINOS_COMPONENTS} SEACASExodus)
  endif()
  find_package(Trilinos REQUIRED COMPONENTS ${NIMBLE_TRILINOS_COMPONENTS})
  target_link_libraries(nimble PUBLIC Trilinos::Trilinos)
else()
  message(STATUS "Trilinos is NOT enabled.")
  set(NIMBLE_HAVE_TRILINOS FALSE)
endif()

# Optional Kokkos dependency
# If trilinos is used, it comes with Kokkos so we don't need standalone
if(NOT NimbleSM_ENABLE_TRILINOS AND NimbleSM_ENABLE_KOKKOS)
  if (Kokkos_ROOT)
    #
    # Search paths specified in the <PackageName>_ROOT
    #
    # See https://cmake.org/cmake/help/latest/command/find_package.html#search-procedure
    #
    find_package(Kokkos REQUIRED)
    set(NIMBLE_HAVE_KOKKOS TRUE)
    set(KOKKOS_HEADER_DIR ${Kokkos_ROOT}/include)
    target_link_libraries(nimble PUBLIC Kokkos::kokkos)
  elseif (KOKKOS_ROOT_DIR)
    message("-- Compiling with Kokkos")
    add_subdirectory(${KOKKOS_ROOT_DIR} ./kokkos)
    include_directories(${Kokkos_INCLUDE_DIRS_RET})
    set(NIMBLE_KOKKOS_LIBS kokkos)
    target_link_libraries(nimble PUBLIC Kokkos::kokkoscore)
    set(NIMBLE_HAVE_KOKKOS TRUE)
  else()
    find_package(Kokkos REQUIRED)
    set(NIMBLE_KOKKOS_LIBS ${Kokkos_LIBRARIES} ${Kokkos_TPL_LIBRARIES})
    set(NIMBLE_HAVE_KOKKOS TRUE)
    message(STATUS "Compiling with Kokkos:")
    message("       Kokkos_DIR: ${Kokkos_DIR}")
    message("       Kokkos_LIBRARIES: ${Kokkos_LIBRARIES}")
    message("       Kokkos_TPL_LIBRARIES: ${Kokkos_TPL_LIBRARIES}")
    message("       Kokkos_CXX_FLAGS: ${Kokkos_CXX_FLAGS}")
    target_link_libraries(nimble PUBLIC Kokkos::kokkos)
  endif()
  target_compile_definitions(nimble PUBLIC -DNIMBLE_HAVE_KOKKOS)
endif()

if(NOT NimbleSM_ENABLE_TRILINOS AND NOT NimbleSM_ENABLE_KOKKOS)
  message("-- Kokkos is NOT enabled.")
  set(NIMBLE_HAVE_KOKKOS FALSE)
  set(NIMBLE_KOKKOS_LIBS)
endif()

# Exodus library dependencies
message(STATUS "HAVE_TRILINOS: ${HAVE_TRILINOS}")
if(NOT NimbleSM_ENABLE_TRILINOS OR NOT NimbleSM_USE_TRILINOS_EXODUS)
  find_package(Exodus REQUIRED)
  set(NIMBLE_HAVE_EXODUS TRUE)
  target_compile_definitions(nimble PRIVATE NIMBLE_HAVE_EXODUS)
  target_link_libraries(nimble PRIVATE Exodus::ExodusII)
  if(CRAY)
    target_link_libraries(nimble PRIVATE hdf5_hl hdf5 z)
  endif()
endif()

# Optional ArborX dependency
if (NimbleSM_ENABLE_ARBORX)
  message(STATUS "Compiling with ArborX:")
  find_package(ArborX REQUIRED)
  #
  add_definitions(-DNIMBLE_HAVE_ARBORX)
  set(NIMBLE_HAVE_ARBORX TRUE)
  target_link_libraries(nimble PRIVATE ArborX::ArborX)
  target_compile_definitions(nimble PUBLIC NIMBLE_HAVE_ARBORX)
  #
  message("       ArborX_DIR: ${ArborX_DIR}")
  if (NIMBLE_HAVE_MPI)
    message("       Using parallel collision detection (NIMBLE_HAVE_MPI=ON)")
  else()
    message("       Using only serial collision detection (NIMBLE_HAVE_MPI=OFF)")
  endif()
else()
  message(STATUS "ArborX is NOT enabled.")
  SET(NIMBLE_HAVE_ARBORX FALSE)
  set(ArborX_LIBS)
endif()

include_directories(src)

add_subdirectory(tpl)
add_subdirectory(src)
add_subdirectory(test)

if(NimbleSM_ENABLE_UNIT_TESTS)
  add_subdirectory(unit_tests)
endif()

install(TARGETS nimble EXPORT NimbleSMTargets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(EXPORT NimbleSMTargets NAMESPACE NimbleSM:: DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/NimbleSM
  EXPORT_LINK_INTERFACE_LIBRARIES)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/NimbleSMSettings.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/NimbleSMSettings.cmake
  @ONLY)

include(CMakePackageConfigHelpers)
configure_package_config_file(cmake/NimbleSMConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/NimbleSMConfig.cmake
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/NimbleSM)

install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/NimbleSMConfig.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/NimbleSMSettings.cmake
  ${CMAKE_SOURCE_DIR}/cmake/FindDARMA.cmake
  ${CMAKE_SOURCE_DIR}/cmake/FindExodus.cmake
  ${CMAKE_SOURCE_DIR}/cmake/FindQThreads.cmake
  ${CMAKE_SOURCE_DIR}/cmake/FindTrilinos.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/NimbleSM)
