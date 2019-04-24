target_sources(pcm-object
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/tsless_cavity.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/tsless_lapack.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/tsless_precision.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/tsless_symmetry.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/tsless_weight_function.f90
  )

