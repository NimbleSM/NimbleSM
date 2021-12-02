
find_package(Trilinos CONFIG COMPONENTS ${${CMAKE_FIND_PACKAGE_NAME}_FIND_COMPONENTS})


# MESSAGE("-- Compiling with Trilinos:")
# MESSAGE("       Trilinos_DIR: ${Trilinos_DIR}")
# MESSAGE("       Trilinos_VERSION: ${Trilinos_VERSION}")
# MESSAGE("       Trilinos_TPL_LIST: ${Trilinos_TPL_LIST}")
# MESSAGE("       Trilinos_PACKAGE_LIST: ${Trilinos_PACKAGE_LIST}")
# message(STATUS "Components: ${Trilinos_FIND_COMPONENTS}")
# MESSAGE("       Trilinos_TPL_INCLUDE_DIRS: ${Trilinos_TPL_INCLUDE_DIRS}")
# MESSAGE("       Trilinos_TPL_LIBRARIES: ${Trilinos_TPL_LIBRARIES}")
# MESSAGE("       Trilinos_TPL_LIBRARY_DIRS: ${Trilinos_TPL_LIBRARY_DIRS}")
# MESSAGE("       Trilinos_BUILD_SHARED_LIBS: ${Trilinos_BUILD_SHARED_LIBS}")
# MESSAGE("       Trilinos_CXX_COMPILER_FLAGS: ${Trilinos_CXX_COMPILER_FLAGS}")
# MESSAGE("       Trilinos_LD_FLAGS: ${Trilinos_LD_FLAGS}")
# MESSAGE("       Trilinos_LIBRARIES: ${Trilinos_LIBRARIES}")
# MESSAGE("       Trilinos_INCLUDE_DIRS: ${Trilinos_INCLUDE_DIRS}")
# MESSAGE("       Trilinos_LIBRARY_DIRS: ${Trilinos_LIBRARY_DIRS}")
string(REPLACE " " ";" _trilinos_flags ${Trilinos_CXX_COMPILER_FLAGS})
if (${Trilinos_LD_FLAGS})
  string(REPLACE " " ";" _trilinos_link_flags ${Trilinos_LD_FLAGS})
endif()

# Trilinos doesn't add -arch=sm_x flags to link flags, so regex them and add manually
foreach(_flag IN LISTS _trilinos_flags)
  if (_flag MATCHES "^--?arch=.+")
    list(APPEND _trilinos_link_flags "${_flag}")
  elseif (_flag MATCHES "^--fopenmp") # fopenmp doens't get propagated either -- this might fail if the openmp flag is different
    list(APPEND _trilinos_link_flags "${_flag}")
  endif()
endforeach()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Trilinos
      REQUIRED_VARS Trilinos_DIR Trilinos_LIBRARIES Trilinos_INCLUDE_DIRS
    )

# Trilinos doesn't use imported targets so add our own here
# Defined as an interface target because there is no explicit library "Trilinos"
# it is a collection of packages
add_library(Trilinos::Trilinos INTERFACE IMPORTED)

set_target_properties(Trilinos::Trilinos PROPERTIES
    INTERFACE_LINK_LIBRARIES "${Trilinos_LIBRARIES};${Trilinos_TPL_LIBRARIES}"
    INTERFACE_INCLUDE_DIRECTORIES "${Trilinos_INCLUDE_DIRS};${Trilinos_TPL_INCLUDE_DIRS}"
    INTERFACE_COMPILE_OPTIONS "${_trilinos_flags}"
    INTERFACE_LINK_OPTIONS "${_trilinos_link_flags}"
    )

set_property(TARGET Trilinos::Trilinos APPEND PROPERTY INTERFACE_LINK_OPTIONS "${_trilinos_link_flags}")
