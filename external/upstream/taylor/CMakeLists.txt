find_package(taylor 1.0.0 CONFIG QUIET)

if(TARGET taylor::taylor)
  get_property(_loc TARGET taylor::taylor PROPERTY INTERFACE_INCLUDE_DIRECTORIES)
  message(STATUS "Using taylor: ${_loc} (found version ${taylor_VERSION})")
  add_library(taylor_external INTERFACE)  # dummy
else()
  include(ExternalProject)
  message(STATUS "Suitable taylor could not be located. Fetching and building!")
  ExternalProject_Add(taylor_external
    URL
      https://github.com/uekstrom/libtaylor/archive/master.zip
    CMAKE_ARGS
      -DCMAKE_INSTALL_PREFIX=${STAGED_INSTALL_PREFIX}
      -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
      -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
      -DCMAKE_CXX_STANDARD=${CMAKE_CXX_STANDARD}
    LOG_DOWNLOAD
      1
    LOG_UPDATE
      1
    LOG_CONFIGURE
      1
    LOG_BUILD
      1
    LOG_INSTALL
      1
    LOG_TEST
      1
    )

  set(taylor_DIR ${STAGED_INSTALL_PREFIX}/share/cmake/taylor CACHE PATH "path to internally built taylorConfig.cmake" FORCE)
endif()
