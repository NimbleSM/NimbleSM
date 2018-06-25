
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/multiscale_anisotropic_brick.in
               ${CMAKE_CURRENT_BINARY_DIR}/multiscale_anisotropic_brick.in COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/multiscale_anisotropic_brick.g
               ${CMAKE_CURRENT_BINARY_DIR}/multiscale_anisotropic_brick.g COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/microstructure_with_fibers_rve.g
               ${CMAKE_CURRENT_BINARY_DIR}/microstructure_with_fibers_rve.g COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/post_process.py
               ${CMAKE_CURRENT_BINARY_DIR}/post_process.py COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/force_displacement.gnuplot
               ${CMAKE_CURRENT_BINARY_DIR}/force_displacement.gnuplot COPYONLY)

