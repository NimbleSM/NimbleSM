#
# This file is only processed when NIMBLE_HAVE_UQ is set to ON
#
set(NIMBLE_PUBLIC_HEADERS ${NIMBLE_PUBLIC_HEADERS}
      ${CMAKE_CURRENT_LIST_DIR}/nimble_uq_block.h
      ${CMAKE_CURRENT_LIST_DIR}/nimble_uq_model_data.h
      ${CMAKE_CURRENT_LIST_DIR}/nimble_uq.h )
#
target_sources(nimble PRIVATE
 ${CMAKE_CURRENT_LIST_DIR}/nimble_uq_block.cc
 ${CMAKE_CURRENT_LIST_DIR}/nimble_uq_model_data.cc
 ${CMAKE_CURRENT_LIST_DIR}/nimble_uq.cc
)

