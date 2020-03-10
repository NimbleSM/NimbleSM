
find_library(QThreads_LIBRARY qthread PATH ${QThreads_DIR} $ENV{QThreads_DIR} PATH_SUFFIXES lib)

find_path(QThreads_INCLUDE_DIRECTORY qthread.h PATH ${QThreads_DIR} $ENV{QThreads_DIR} PATH_SUFFIXES include)

mark_as_advanced(QThreads_LIBRARY)
mark_as_advanced(QThreads_INCLUDE_DIRECTORY)

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(QThreads
    REQUIRED_VARS
    QThreads_LIBRARY
    QThreads_INCLUDE_DIRECTORY
    )

add_library(QThreads::QThreads UNKNOWN IMPORTED)

set_target_properties(QThreads::QThreads PROPERTIES
    IMPORTED_LOCATION "${QThreads_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${QThreads_INCLUDE_DIRECTORY}")