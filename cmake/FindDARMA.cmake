#.rst
# DARMA
# -----
#
# Searches for an installation of the DARMA tasking library.
# Sets the following variables:
# DARMA_BACKEND_INCLUDE_DIRS
# DARMA_BACKEND_LIBRARY_DIRS
# DARMA_BACKEND_CXX_COMPILER
# DARMA_BACKEND_CXX_FLAGS
# DARMA_BACKEND_LIBRARIES
# DARMA_INCLUDE_DIRS
# DARMA_LIBRARIES
#
# Specifying the following variables can impact the behavior of FindDarma:
# DISABLE_SPACK - Don't use the spack package manager for finding DARMA
# DARMA_INCLUDE_DIRS - FindDarma will first look in here for darma.h rather than using spack
# DARMA_LIBRARY_DIRS - FindDarma will first look in here for the DARMA library files rather than using spack or Config.cmake
# DARMA_FRONTEND_INCLUDE_DIR - FindDarma will use this for the DARMA frontend include dir rather than trying to find it
# DARMA_BACKEND_ROOT - location to search for the DARMA .cmake configure scripts rather than using spack
#

# Convert a component name to a config file to search for. Handles special cases such as Threads being SimpleBackend
MACRO(CONFIG_FROM_COMPONENT COMPONENT)
  if(${COMPONENT} STREQUAL "THREADS")
    SET(_DARMA_CONFIG "DarmaSimpleBackend")
  else()
    SET(_DARMA_CONFIG "Darma${COMPONENT}Backend")
  endif()
ENDMACRO()

MACRO(GET_COMPILER_PATH SPACK_COMPILER_SPEC)
  EXECUTE_PROCESS(COMMAND ${_DARMA_SPACK_PROGRAM} compiler info ${SPACK_COMPILER_SPEC}
                  COMMAND awk "\$1 == \"cxx\" {print \$3}"
                  OUTPUT_VARIABLE _DARMA_COMPILER
                  OUTPUT_STRIP_TRAILING_WHITESPACE)
  MESSAGE("Compiler path: ${_DARMA_COMPILER}")

ENDMACRO()

# Get the correct compiler string
MACRO(GET_COMPILER_VERSION DARMA_SPACK_PACKAGE)
  SET(_AWK_SPLIT "%")
  EXECUTE_PROCESS(COMMAND ${_DARMA_SPACK_PROGRAM} find --show-full-compiler -f ${DARMA_SPACK_PACKAGE}
                  COMMAND awk "$1 ~ /^${DARMA_SPACK_PACKAGE}/ {print \$1}"
                  COMMAND awk -F${_AWK_SPLIT} "{print \$2}"
                  OUTPUT_VARIABLE _DARMA_SPACK_COMPILER_SPEC
                  OUTPUT_STRIP_TRAILING_WHITESPACE)
  MESSAGE("Found spack compiler specs: ${_DARMA_SPACK_COMPILER_SPEC}")
ENDMACRO()

# Only valid for spack
MACRO(FIND_DEPENDENT_FRONTEND COMPONENT)
  STRING(TOLOWER ${COMPONENT} _COMPONENT_SUBNAME)
  SET(_DARMA_${COMPONENT}_NAME "darma-${_COMPONENT_SUBNAME}")

  if (NOT (_DARMA_SPACK_PROGRAM-NOTFOUND OR DARMA_DISABLE_SPACK))
    EXECUTE_PROCESS(COMMAND ${_DARMA_SPACK_PROGRAM} dependencies -i ${_DARMA_${COMPONENT}_NAME}
                    COMMAND awk "(NR > 2) && (/darma-frontend/) { print \$2 }"
                    OUTPUT_VARIABLE _DARMA_SPACK_${COMPONENT}_FRONTEND
                    OUTPUT_STRIP_TRAILING_WHITESPACE)
    MESSAGE("Found dependent frontend for darma component ${COMPONENT} package: ${_DARMA_SPACK_${COMPONENT}_FRONTEND}")

    EXECUTE_PROCESS(COMMAND ${_DARMA_SPACK_PROGRAM} location -i ${_DARMA_SPACK_${COMPONENT}_FRONTEND}
        OUTPUT_VARIABLE _DARMA_SPACK_FRONTEND_DIR
        OUTPUT_STRIP_TRAILING_WHITESPACE)
  endif()
ENDMACRO()

#Set backend vars
MACRO(FIND_BACKEND COMPONENT)
  LIST(APPEND _DARMA_REQUIRED_VARS DARMA_${COMPONENT}_LIBRARY)
  STRING(TOLOWER ${COMPONENT} _COMPONENT_SUBNAME)
  SET(_DARMA_${COMPONENT}_NAME "darma-${_COMPONENT_SUBNAME}")

  if (NOT (_DARMA_SPACK_PROGRAM-NOTFOUND OR DARMA_DISABLE_SPACK))
    EXECUTE_PROCESS(COMMAND ${_DARMA_SPACK_PROGRAM} location -i ${_DARMA_${COMPONENT}_NAME}
        OUTPUT_VARIABLE _DARMA_SPACK_${COMPONENT}_DIR
        OUTPUT_STRIP_TRAILING_WHITESPACE)
  endif()

  # Get the cmake files for this
  # Sets:
  # DARMA_BACKEND_INCLUDE_DIRS
  # DARMA_BACKEND_LIBRARY_DIRS
  # DARMA_BACKEND_CXX_COMPILER
  # DARMA_BACKEND_CXX_FLAGS
  # DARMA_BACKEND_LIBRARIES
  CONFIG_FROM_COMPONENT(${COMPONENT})
  MESSAGE("Using config ${_DARMA_CONFIG}")
  FIND_PACKAGE(${_DARMA_CONFIG} PATHS ${DARMA_BACKEND_ROOT} ${_DARMA_SPACK_${COMPONENT}_DIR})

  if (NOT ${_DARMA_CONFIG}_FOUND)
    MESSAGE(FATAL_ERROR "Backend at ${_DARMA_SPACK_${COMPONENT}_DIR} is not a valid backend, check your install")
  endif()


  SET(DARMA_${COMPONENT}_LIBRARY_DIR ${DARMA_BACKEND_LIBRARY_DIRS})
  LIST(APPEND DARMA_BACKEND_LIBRARY_DIR ${DARMA_${COMPONENT}_LIBRARY_DIR})
  MESSAGE("Darma library dirs: ${DARMA_BACKEND_LIBRARY_DIRS}")

  foreach(LIBRARY ${DARMA_BACKEND_LIBRARIES})
    MESSAGE("Finding library ${LIBRARY}")
    STRING(TOUPPER ${LIBRARY} _LIBRARY)
    FIND_LIBRARY(${_LIBRARY}_LIBRARY
        NAMES ${LIBRARY}
        PATHS ${DARMA_LIBRARY_DIRS} ${DARMA_BACKEND_LIBRARY_DIRS})
    if (${${_LIBRARY}_LIBRARY} STREQUAL "${_LIBRARY}_LIBRARY-NOTFOUND")
      MESSAGE(FATAL_ERROR "Could not find library ${LIBRARY}")
    else()
      MESSAGE("Found library ${${_LIBRARY}_LIBRARY}")
    endif()
    LIST(APPEND DARMA_BACKEND_LIBRARY ${${_LIBRARY}_LIBRARY})


    # Add the main library too
    # TODO: remove this when darma config is fixed
    FIND_LIBRARY(${_LIBRARY}_LIBRARY_MAIN
        NAMES ${LIBRARY}_main
        PATHS ${DARMA_LIBRARY_DIRS} ${DARMA_BACKEND_LIBRARY_DIRS})
    if (${${_LIBRARY}_LIBRARY_MAIN} STREQUAL "${_LIBRARY}_LIBRARY_MAIN-NOTFOUND")
      MESSAGE(WARNING "Could not find library ${LIBRARY}_main. Depending on how DARMA is configured this may not be a problem.")
    else()
      MESSAGE("Found library ${${_LIBRARY}_LIBRARY_MAIN}")
      LIST(APPEND DARMA_BACKEND_LIBRARY_MAIN ${${_LIBRARY}_LIBRARY_MAIN})
    endif()
  endforeach()

  # Try and find the darma_types header to locate backend header directory
  FIND_PATH(DARMA_BACKEND_INCLUDE_DIR
      NAMES darma_types.h
      PATHS  ${DARMA_INCLUDE_DIRS} ${_DARMA_SPACK_${COMPONENT}_DIR}/include ENV ${CPATH}
      DOC "Path to the include directory for the ${COMPONENT} backend.")

  MESSAGE("Backend path: ${DARMA_BACKEND_INCLUDE_DIR}")

  GET_COMPILER_VERSION(${_DARMA_${COMPONENT}_NAME})
  GET_COMPILER_PATH(${_DARMA_SPACK_COMPILER_SPEC})
ENDMACRO()


if (NOT DARMA_DISABLE_SPACK)
  # Append spack location search paths
  FIND_PROGRAM(_DARMA_SPACK_PROGRAM spack)
  SET(DARMA_FRONTEND_NAME "darma-frontend")
  if (NOT _DARMA_SPACK_PROGRAM-NOTFOUND)
    MESSAGE("spack found in ${_DARMA_SPACK_PROGRAM}")
  endif()
endif()

foreach(COMPONENT ${DARMA_FIND_COMPONENTS})
  STRING(TOUPPER ${COMPONENT} _COMPONENT)
  SET(DARMA_USE_${_COMPONENT} 1)


  FIND_DEPENDENT_FRONTEND(${_COMPONENT})
  FIND_BACKEND(${_COMPONENT})
endforeach()

# Deal with multiple backends, can only have one
LIST(LENGTH DARMA_BACKEND_LIBRARY _DARMA_NUM_BACKENDS)
if (NOT _DARMA_NUM_BACKENDS EQUAL 1)
  if (_DARMA_NUM_BACKENDS GREATER 1)
    MESSAGE(FATAL_ERROR "Multiple backends defined: ${DARMA_BACKEND_LIBRARY}")
  else()
    MESSAGE(FATAL_ERROR "No valid backends")
  endif()
endif()

# main darma header always required
SET(_DARMA_REQUIRED_VARS DARMA_FRONTEND_INCLUDE_DIR)
FIND_PATH(DARMA_FRONTEND_INCLUDE_DIR
    NAMES darma.h
    PATHS ${DARMA_INCLUDE_DIRS} ${_DARMA_SPACK_FRONTEND_DIR}/include ENV ${CPATH})


include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(DARMA
    FOUND_VAR DARMA_FOUND
    REQUIRED_VARS _DARMA_REQUIRED_VARS
    )

if (DARMA_FOUND)
  SET(DARMA_INCLUDE_DIRS ${DARMA_FRONTEND_INCLUDE_DIR} ${DARMA_BACKEND_INCLUDE_DIR})
  SET(DARMA_BACKEND_INCLUDE_DIRS ${DARMA_BACKEND_INCLUDE_DIR})
  SET(DARMA_BACKEND_LIBRARY_DIRS ${DARMA_BACKEND_LIBRARY_DIR})
  SET(DARMA_LIBRARIES ${DARMA_BACKEND_LIBRARY} ${DARMA_BACKEND_LIBRARY_MAIN})
  SET(DARMA_BACKEND_LIBRARIES ${DARMA_BACKEND_LIBRARY} ${DARMA_BACKEND_LIBRARY_MAIN})
  SET(DARMA_BACKEND_CXX_COMPILER ${_DARMA_COMPILER})
endif()
