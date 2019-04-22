cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(test-PCMSolver LANGUAGES CXX)

enable_testing()

include(CMakeParseArguments)
macro(add_Catch_test)
  set(oneValueArgs NAME COST)
  set(multiValueArgs LABELS DEPENDS REFERENCE_FILES)
  cmake_parse_arguments(add_Catch_test
    "${options}"
    "${oneValueArgs}"
    "${multiValueArgs}"
    ${ARGN}
    )

  add_test(
    NAME
      ${add_Catch_test_NAME}
    COMMAND
      ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/unit_tests
      [${add_Catch_test_NAME}] --success --out
      ${PROJECT_BINARY_DIR}/${add_Catch_test_NAME}.log --durations yes
    WORKING_DIRECTORY
      ${CMAKE_CURRENT_BINARY_DIR}
    )

  set_tests_properties(${add_Catch_test_NAME}
    PROPERTIES
      LABELS "${add_Catch_test_LABELS}"
    )

  if(add_Catch_test_COST)
    set_tests_properties(${add_Catch_test_NAME}
      PROPERTIES
        COST ${add_Catch_test_COST}
      )
  endif()

  if(add_Catch_test_DEPENDS)
    set_tests_properties(${add_Catch_test_NAME}
      PROPERTIES
        DEPENDS ${add_Catch_test_DEPENDS}
      )
  endif()

  if(add_Catch_test_REFERENCE_FILES)
    file(
      COPY
        ${add_Catch_test_REFERENCE_FILES}
      DESTINATION
        ${CMAKE_CURRENT_BINARY_DIR}
      )
  endif()
endmacro()

find_package(Eigen3 3.3 CONFIG REQUIRED)
message(STATUS "Using Eigen3: ${EIGEN3_ROOT_DIR} (version ${Eigen3_VERSION})")

find_package(PCMSolver CONFIG REQUIRED)
get_property(_loc TARGET PCMSolver::pcm PROPERTY INTERFACE_INCLUDE_DIRECTORIES)
message(STATUS "Using PCMSolver: ${_loc} (found version ${PCMSolver_VERSION})")

set(CMAKE_INCLUDE_CURRENT_DIR ON)

add_executable(unit_tests unit_tests.cpp)

add_subdirectory(input)
add_subdirectory(numerical_quadrature)
add_subdirectory(gepol)
add_subdirectory(dielectric_profile)
add_subdirectory(green)
add_subdirectory(bi_operators)
add_subdirectory(cpcm)
add_subdirectory(iefpcm)
add_subdirectory(utils)
add_subdirectory(benchmark)

target_link_libraries(unit_tests
  PRIVATE
    Eigen3::Eigen
    PCMSolver::pcm
  )

#add_subdirectory(C_host)
#if(TEST_Fortran_API)
#  add_subdirectory(Fortran_host)
#endif()
#if(BUILD_STANDALONE)
#  add_subdirectory(standalone)
#endif()
