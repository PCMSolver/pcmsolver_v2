find_package(Eigen3 3.3 CONFIG QUIET)

if(TARGET Eigen3::Eigen)
  message(STATUS "Using Eigen3: ${EIGEN3_ROOT_DIR} (version ${Eigen3_VERSION})")
  add_library(Eigen3_external INTERFACE)  # dummy
else()
  include(ExternalProject)
  message(STATUS "Suitable Eigen3 could not be located. Fetching and building!")
  ExternalProject_Add(Eigen3_external
    URL
      https://github.com/eigenteam/eigen-git-mirror/archive/3.3.7.zip
    CMAKE_ARGS
      -DCMAKE_INSTALL_PREFIX=${STAGED_INSTALL_PREFIX}
      -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
      -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER} # Unused but probed
      -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
      -DCMAKE_CXX_STANDARD=${CMAKE_CXX_STANDARD}
      -DBUILD_TESTING=OFF
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

  set(Eigen3_DIR ${STAGED_INSTALL_PREFIX}/share/eigen3/cmake CACHE PATH "path to internally built Eigen3Config.cmake" FORCE)
endif()
