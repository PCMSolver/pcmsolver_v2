cmake_minimum_required(VERSION 3.5 FATAL_ERROR)

project(test-PCMSolver LANGUAGES CXX C Fortran)

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
get_property(_loc TARGET PCMSolver::pcm PROPERTY LOCATION)
message(STATUS "Using PCMSolver: ${_loc} (found version ${PCMSolver_VERSION})")

find_package(PythonInterp REQUIRED)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

add_executable(unit_tests unit_tests.cpp)

# To make tests accessible from the parent project:
#   1. We save the testing info for each folder in the parent binary directory.
#   2. We install CTestTestfile.cmake in the parent binary directory.

add_subdirectory(input ${PARENT_BINARY_DIR}/tests/input)
add_subdirectory(numerical_quadrature ${PARENT_BINARY_DIR}/tests/numerical_quadrature)
add_subdirectory(gepol ${PARENT_BINARY_DIR}/tests/gepol)
add_subdirectory(dielectric_profile ${PARENT_BINARY_DIR}/tests/dielectric_profile)
add_subdirectory(green ${PARENT_BINARY_DIR}/tests/green)
add_subdirectory(bi_operators ${PARENT_BINARY_DIR}/tests/bi_operator)
add_subdirectory(cpcm ${PARENT_BINARY_DIR}/tests/cpcm)
add_subdirectory(iefpcm ${PARENT_BINARY_DIR}/tests/iefpcm)
add_subdirectory(utils ${PARENT_BINARY_DIR}/tests/utils)
# Not sure what to do about these, actually
#add_subdirectory(benchmark ${PARENT_BINARY_DIR}/tests/benchmark)

# Little hack to avoid prefixing every include in every source file with "PCMSolver"
target_include_directories(unit_tests
  SYSTEM PRIVATE
    ${PCMSolver_INCLUDE_DIRS}/PCMSolver
  )

target_link_libraries(unit_tests
  PRIVATE
    Eigen3::Eigen
    PCMSolver::pcm
  )

install(
  FILES
    ${CMAKE_CURRENT_BINARY_DIR}/CTestTestfile.cmake
  DESTINATION
    ${PARENT_BINARY_DIR}
  )

# API tests
add_subdirectory(C_host ${PARENT_BINARY_DIR}/tests/C_host)
add_subdirectory(Fortran_host ${PARENT_BINARY_DIR}/tests/Fortran_host)
add_subdirectory(standalone ${PARENT_BINARY_DIR}/tests/standalone)
