#
# CMakeLists.txt file for NimbleSM
#

message("CMakeLists.txt - NimbleSM")

cmake_minimum_required(VERSION 3.12)
project(NimbleSM VERSION 0.1.0 LANGUAGES CXX)

include(GNUInstallDirs)

if (NIMBLE_NVIDIA_BUILD)
  add_definitions(-DNIMBLE_NVIDIA_BUILD)
endif()

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

enable_testing()

set(SRC_DIR ${CMAKE_SOURCE_DIR}/src)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

option(HAVE_UQ "Whether to use enamble UQ sampling" OFF)
option(HAVE_BVH "Whether to use BVH for parallel asynchronous contact" OFF)
option(HAVE_TRILINOS "Whether to use Trilinos" OFF)
option(HAVE_KOKKOS "Whether to use Kokkos for parallel execution" OFF)
option(USE_PURE_MPI "Whether to use Pure MPI (no Trilinos)" OFF)
option(HAVE_QTHREADS "Whether to use QThreads" OFF)
option(NIMBLE_NVIDIA_BUILD "Whether to build Nimble with CUDA" OFF)
option(NIMBLE_ENABLE_UNIT_TESTS "Whether to build Nimble with GTest testing" OFF)
option(TIME_CONTACT "Whether to time contact" OFF)
option(HAVE_ARBORX "Whether to use ArborX" OFF)

add_library(nimble)
add_library(nimble::nimble ALIAS nimble)

target_include_directories(nimble PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

# Optional functionality for UQ
if (HAVE_UQ)
  set(NIMBLE_HAVE_UQ TRUE)
  message(STATUS "Compiling with UQ")
  target_compile_definitions(nimble PUBLIC NIMBLE_HAVE_UQ)
endif()

# Optional functionality for contact
if (TIME_CONTACT)
  set(NIMBLE_TIME_CONTACT TRUE)
  message(STATUS "Compiling with contact timing")
  target_compile_definitions(nimble PUBLIC NIMBLE_TIME_CONTACT)
endif()

# Optional dependency for pure MPI build (no Trilinos)
if (USE_PURE_MPI)
  set(NIMBLE_HAVE_MPI TRUE)
    message(STATUS "Compiling with MPI")
  find_package(MPI REQUIRED)
  target_link_libraries(nimble PUBLIC MPI::MPI_CXX)
  target_compile_definitions(nimble PUBLIC NIMBLE_HAVE_MPI)
endif()

# Optional BVH dependency
if (HAVE_BVH)
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
if(HAVE_TRILINOS)
  option(USE_TRILINOS_EXODUS "Whether to use Exodus from Trilinos" ON)
  # Trilinos brings in Trilinos, Kokkos, and (optionally) Exodus
  target_compile_definitions(nimble PUBLIC -DNIMBLE_HAVE_TRILINOS)
  set(NIMBLE_HAVE_TRILINOS TRUE)
  target_compile_definitions(nimble PUBLIC -DNIMBLE_HAVE_KOKKOS)
  set(NIMBLE_HAVE_KOKKOS TRUE)
  if(USE_TRILINOS_EXODUS)
    target_compile_definitions(nimble PUBLIC -DNIMBLE_HAVE_EXODUS)
    set(NIMBLE_HAVE_EXODUS TRUE)
  endif()
  message(STATUS "Trilinos is enabled.")
  set(Trilinos_DIR CACHE PATH "Path to Trilinos install directory")
  set(NIMBLE_TRILINOS_COMPONENTS KokkosCore TpetraCore)
  if(USE_TRILINOS_EXODUS)
    set(NIMBLE_TRILINOS_COMPONENTS ${NIMBLE_TRILINOS_COMPONENTS} SEACASExodus)
  endif()
  find_package(Trilinos REQUIRED COMPONENTS ${NIMBLE_TRILINOS_COMPONENTS})
  target_link_libraries(nimble PUBLIC Trilinos::Trilinos)
else()
  message(STATUS "Trilinos is NOT enabled.")
  set(NIMBLE_HAVE_TRILINOS FALSE)
endif()

# Qthreads library dependencies
if (HAVE_QTHREADS)
  find_package(QThreads REQUIRED)
  target_link_libraries(nimble PUBLIC QThreads::QThreads)
  target_compile_definitions(nimble PUBLIC -DNIMBLE_HAVE_QTHREADS)
  set(NIMBLE_HAVE_QTHREADS TRUE)
endif()

# Optional Kokkos dependency
# If trilinos is used, it comes with Kokkos so we don't need standalone
if(NOT HAVE_TRILINOS AND HAVE_KOKKOS)
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

if(NOT HAVE_TRILINOS AND NOT HAVE_KOKKOS)
  message("-- Kokkos is NOT enabled.")
  set(NIMBLE_HAVE_KOKKOS FALSE)
  set(NIMBLE_KOKKOS_LIBS)
endif()

# Exodus library dependencies
message(STATUS "HAVE_TRILINOS: ${HAVE_TRILINOS}")
if(NOT HAVE_TRILINOS OR NOT USE_TRILINOS_EXODUS)
  find_package(Exodus REQUIRED)
  set(NIMBLE_HAVE_EXODUS TRUE)
  target_compile_definitions(nimble PRIVATE NIMBLE_HAVE_EXODUS)
  target_link_libraries(nimble PRIVATE Exodus::ExodusII)
  if(CRAY)
    target_link_libraries(nimble PRIVATE hdf5_hl hdf5 z)
  endif()
endif()

# Optional ArborX dependency
if (HAVE_ARBORX)
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

add_subdirectory(src)
add_subdirectory(test)

if(NIMBLE_ENABLE_UNIT_TESTS)
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
