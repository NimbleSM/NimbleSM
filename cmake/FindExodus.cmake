
# First find an Exodus that has a cmake config
find_package(SEACASExodus CONFIG QUIET)
if (SEACASExodus_FOUND)
  find_library(exodus_LIBRARY NAMES exodus PATHS "${SEACASExodus_LIBRARY_DIRS}")
  mark_as_advanced(exodus_LIBRARY)

  find_path(exodus_INCLUDE_DIR NAMES exodusII.h PATHS "${SEACASExodus_INCLUDE_DIRS}")
  mark_as_advanced(exodus_INCLUDE_DIR)
endif()


# Version info
if (exodus_INCLUDE_DIR)
  file(STRINGS "${exodus_INCLUDE_DIR}/exodusII.h" _exodus_version_lines REGEX "#define[ \t]+EX_API_VERS")
  if("${_exodus_version_lines}" MATCHES "#define[ \t]+EX_API_VERS[ \t]+\([0-9]\)+\.\([0-9]+\)f")
    set(Exodus_VERSION_MAJOR "${CMAKE_MATCH_1}")
    set(Exodus_VERSION_MINOR "${CMAKE_MATCH_2}")
    set(Exodus_VERSION "${Exodus_VERSION_MAJOR}.${Exodus_VERSION_MINOR}")
    set(Exodus_VERSION_STRING "${Exodus_VERSION}")
  endif()
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Exodus REQUIRED_VARS exodus_LIBRARY exodus_INCLUDE_DIR VERSION_VAR Exodus_VERSION)

message(STATUS "exodus libraries: ${SEACASExodus_TPL_LIBRARIES}")

if (Exodus_FOUND)
  set(Exodus_INCLUDE_DIRS "${exodus_INCLUDE_DIR}")
  set(Exodus_LIBRARIES "${exodus_LIBRARY}")

  if (NOT TARGET Exodus::ExodusII)
    add_library(Exodus::ExodusII UNKNOWN IMPORTED)
    set_target_properties(Exodus::ExodusII PROPERTIES
      IMPORTED_LOCATION "${exodus_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${exodus_INCLUDE_DIR}"#;${SEACASExodus_TPL_INCLUDE_DIRS}"
      #INTERFACE_LINK_LIBRARIES "${SEACASExodus_TPL_LIBRARIES}"
    )
  endif()
endif()

