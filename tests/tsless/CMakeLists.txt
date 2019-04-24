target_sources(unit_tests
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/tsless_HF.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/tsless_NH3.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/tsless_point.cpp
  )

# tsless_HF.cpp test
add_Catch_test(
  NAME
    tsless_HF
  LABELS
    tsless
    tsless_HF
  )

# tsless_point.cpp test
add_Catch_test(
  NAME
    tsless_point
  LABELS
    tsless
    tsless_point
  )

# tsless_NH3.cpp test
add_Catch_test(
  NAME
    tsless_NH3
  LABELS
    tsless
    tsless_NH3
  )

