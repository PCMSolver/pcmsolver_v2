cmake_minimum_required(VERSION 3.5 FATAL_ERROR)
project(PCMSolver
  LANGUAGES CXX Fortran
  VERSION 2.0.0
  )

include(GNUInstallDirs)

# Offer the user the choice of overriding the installation directories
set(INSTALL_LIBDIR ${CMAKE_INSTALL_LIBDIR} CACHE PATH "Installation directory for libraries")
set(INSTALL_BINDIR ${CMAKE_INSTALL_BINDIR} CACHE PATH "Installation directory for executables")
set(INSTALL_INCLUDEDIR ${CMAKE_INSTALL_INCLUDEDIR} CACHE PATH "Installation directory for header files")
if(WIN32 AND NOT CYGWIN)
  set(DEF_INSTALL_CMAKEDIR CMake)
else()
  set(DEF_INSTALL_CMAKEDIR share/cmake/${PROJECT_NAME})
endif()
set(INSTALL_CMAKEDIR ${DEF_INSTALL_CMAKEDIR} CACHE PATH "Installation directory for CMake files")

# Report to user
foreach(p LIB BIN INCLUDE CMAKE)
  file(TO_NATIVE_PATH ${CMAKE_INSTALL_PREFIX}/${INSTALL_${p}DIR} _path )
  message(STATUS "Installing ${p} components to ${_path}")
  unset(_path)
endforeach()

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/${INSTALL_LIBDIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/${INSTALL_LIBDIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/${INSTALL_BINDIR})

# <<< Dependencies >>>

find_package(Eigen3 3.3 CONFIG REQUIRED)
message(STATUS "Using Eigen3: ${EIGEN3_ROOT_DIR} (version ${Eigen3_VERSION})")

find_package(taylor 1.0.0 CONFIG REQUIRED)
get_property(_loc TARGET taylor::taylor PROPERTY INTERFACE_INCLUDE_DIRECTORIES)
message(STATUS "Using taylor: ${_loc} (found version ${taylor_VERSION})")

find_package(ZLIB REQUIRED)

find_package(Boost REQUIRED) # REMOVE!!!!

add_library(pcm-object OBJECT "")

target_include_directories(pcm-object
  PRIVATE
    ${PROJECT_SOURCE_DIR}
    ${PROJECT_SOURCE_DIR}/api
    ${PROJECT_BINARY_DIR}
  )

target_link_libraries(pcm-object
  PRIVATE
    Eigen3::Eigen
    taylor::taylor
    Boost::boost # REMOVE!!!
    ZLIB::ZLIB
  )

add_subdirectory(bi_operators)
add_subdirectory(cavity)
add_subdirectory(green)
add_subdirectory(interface)
add_subdirectory(pedra)
add_subdirectory(solver)
add_subdirectory(utils)
add_subdirectory(tsless)

add_library(pcm $<TARGET_OBJECTS:pcm-object>)

include(GenerateExportHeader)
generate_export_header(pcm
  BASE_NAME "PCMSolver"
  EXPORT_MACRO_NAME "PCMSolver_EXPORT"
  EXPORT_FILE_NAME "PCMSolverExport.h"
  DEPRECATED_MACRO_NAME "PCMSolver_DEPRECATED"
  NO_EXPORT_MACRO_NAME "PCMSolver_NO_EXPORT"
  STATIC_DEFINE "PCMSolver_STATIC_DEFINE"
  NO_DEPRECATED_MACRO_NAME "PCMSolver_NO_DEPRECATED"
  DEFINE_NO_DEPRECATED
  )

list(APPEND _public_headers
  ${PROJECT_SOURCE_DIR}/api/pcmsolver.h
  ${PROJECT_SOURCE_DIR}/api/pcmsolver.f90
  ${PROJECT_SOURCE_DIR}/api/PCMInput.h
  ${PROJECT_BINARY_DIR}/PCMSolverExport.h
  )

add_subdirectory(tools)

target_compile_definitions(pcm
  INTERFACE
    $<INSTALL_INTERFACE:USING_PCMSolver>
  )

target_include_directories(pcm
  INTERFACE
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/api>
    $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>
    $<INSTALL_INTERFACE:${INSTALL_INCLUDEDIR}>
  )

target_link_libraries(pcm
  INTERFACE
    ZLIB::ZLIB
  )

set_target_properties(pcm
  PROPERTIES
    SOVERSION ${PROJECT_VERSION_MAJOR}
    OUTPUT_NAME "pcm"
    DEBUG_POSTFIX "_d"
    PUBLIC_HEADER "${_public_headers}"
    MACOSX_RPATH ON
  )

add_subdirectory(bin)

install(
  TARGETS
    pcm
  EXPORT
    "PCMSolverTargets"
  RUNTIME DESTINATION
    ${INSTALL_LIBDIR}
    COMPONENT bin
  LIBRARY DESTINATION
    ${INSTALL_LIBDIR}
    COMPONENT lib
  ARCHIVE DESTINATION
    ${INSTALL_LIBDIR}
    COMPONENT lib
  PUBLIC_HEADER
    DESTINATION ${INSTALL_INCLUDEDIR}/PCMSolver
    COMPONENT dev
  )

# <<<  Export Config  >>>

include(CMakePackageConfigHelpers)

write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion
  )

configure_package_config_file(
  ${PROJECT_SOURCE_DIR}/${PROJECT_NAME}Config.cmake.in
  ${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
  INSTALL_DESTINATION
    ${INSTALL_CMAKEDIR}
  )
install(
  FILES
    ${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
    ${PROJECT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
  DESTINATION
    ${INSTALL_CMAKEDIR}
  )

install(
  EXPORT
    "PCMSolverTargets"
  NAMESPACE
    "PCMSolver::"
  DESTINATION
    ${INSTALL_CMAKEDIR}
  )
