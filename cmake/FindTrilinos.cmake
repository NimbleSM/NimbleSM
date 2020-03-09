
find_package(Trilinos CONFIG COMPONENTS ${${CMAKE_FIND_PACKAGE_NAME}_FIND_COMPONENTS})


# MESSAGE("-- Compiling with Trilinos:")
# MESSAGE("       Trilinos_DIR: ${Trilinos_DIR}")
# MESSAGE("       Trilinos_VERSION: ${Trilinos_VERSION}")
# MESSAGE("       Trilinos_TPL_LIST: ${Trilinos_TPL_LIST}")
# MESSAGE("       Trilinos_PACKAGE_LIST: ${Trilinos_PACKAGE_LIST}")
# message(STATUS "COmponents: ${Trilinos_FIND_COMPONENTS}")
# MESSAGE("       Trilinos_TPL_INCLUDE_DIRS: ${Trilinos_TPL_INCLUDE_DIRS}")
# MESSAGE("       Trilinos_TPL_LIBRARIES: ${Trilinos_TPL_LIBRARIES}")
# MESSAGE("       Trilinos_TPL_LIBRARY_DIRS: ${Trilinos_TPL_LIBRARY_DIRS}")
# MESSAGE("       Trilinos_BUILD_SHARED_LIBS: ${Trilinos_BUILD_SHARED_LIBS}")
# MESSAGE("       Trilinos_CXX_COMPILER_FLAGS: ${Trilinos_CXX_COMPILER_FLAGS}")
#MESSAGE("       Trilinos_LIBRARIES: ${Trilinos_LIBRARIES}")
# MESSAGE("       Trilinos_INCLUDE_DIRS: ${Trilinos_INCLUDE_DIRS}")
# MESSAGE("       Trilinos_LIBRARY_DIRS: ${Trilinos_LIBRARY_DIRS}")


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
    INTERFACE_LINK_OPTIONS "${Trilinos_LD_FLAGS}"
    )