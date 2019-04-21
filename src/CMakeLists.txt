cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(PCMSolver LANGUAGES CXX)

include(GNUInstallDirs)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})

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

# <<< Dependencies >>>

find_package(Eigen3 3.3 CONFIG REQUIRED)
message(STATUS "Using Eigen3: ${EIGEN3_ROOT_DIR} (version ${Eigen3_VERSION})")

find_package(taylor 1.0.0 CONFIG REQUIRED)
get_property(_loc TARGET taylor::taylor PROPERTY INTERFACE_INCLUDE_DIRECTORIES)
message(STATUS "Using taylor: ${_loc} (found version ${taylor_VERSION})")

find_package(ZLIB REQUIRED)

find_package(Boost REQUIRED) # REMOVE!!!!

add_library(pcm "")

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

target_include_directories(pcm
  PRIVATE
    ${PROJECT_SOURCE_DIR}
    ${PROJECT_SOURCE_DIR}/api
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_BINARY_DIR}
  )

target_link_libraries(pcm
  PUBLIC
    Eigen3::Eigen
    taylor::taylor
    Boost::boost # REMOVE!!!
  )

add_subdirectory(cavity)
add_subdirectory(green)
add_subdirectory(interface)
add_subdirectory(bi_operators)
add_subdirectory(pedra)
add_subdirectory(solver)
add_subdirectory(utils)

# <<<  Export Config  >>>

#include(CMakePackageConfigHelpers)
#configure_package_config_file(
#  ${PROJECT_SOURCE_DIR}/${PROJECT_NAME}Config.cmake.in
#  ${PROJECT_BINARY_DIR}/${INSTALL_CMAKEDIR}/${PROJECT_NAME}Config.cmake
#  INSTALL_DESTINATION
#    ${INSTALL_CMAKEDIR}
#  )
#install(
#  FILES
#    ${PROJECT_BINARY_DIR}/${INSTALL_CMAKEDIR}/${PROJECT_NAME}Config.cmake
#    ${PROJECT_BINARY_DIR}/${INSTALL_CMAKEDIR}/${PROJECT_NAME}ConfigVersion.cmake
#  DESTINATION
#    ${INSTALL_CMAKEDIR}
#  )
#
#install(
#  EXPORT
#    "PCMSolverTargets"
#  NAMESPACE
#    "PCMSolver::"
#  DESTINATION
#    ${INSTALL_CMAKEDIR}
#  )
