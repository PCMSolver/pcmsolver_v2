# Executables for testing separately from the unit test framework or for profiling purposes

add_executable(run_pcm run_pcm.cpp)

target_link_libraries(run_pcm
  PRIVATE
    Eigen3::Eigen
    pcm
  )

if(ENABLE_GPERFTOOLS)
  target_link_libraries(run_pcm ${GPERFTOOLS_LIBRARIES})
endif()

# Prepare RPATH

file(RELATIVE_PATH _rel ${CMAKE_INSTALL_PREFIX}/${INSTALL_BINDIR} ${CMAKE_INSTALL_PREFIX})
if(APPLE)
  set(_rpath "@loader_path/${_rel}")
else()
  set(_rpath "\$ORIGIN/${_rel}")
endif()
file(TO_NATIVE_PATH "${_rpath}/${INSTALL_LIBDIR}" _RPATH)

set_target_properties(run_pcm
  PROPERTIES
    MACOSX_RPATH ON
    SKIP_BUILD_RPATH OFF
    BUILD_WITH_INSTALL_RPATH OFF
    INSTALL_RPATH "${_RPATH}"
    INSTALL_RPATH_USE_LINK_PATH ON
)

install(
  TARGETS
    run_pcm
  RUNTIME DESTINATION
    ${INSTALL_BINDIR}
    COMPONENT bin
  )
