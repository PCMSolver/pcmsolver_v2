target_sources(pcm-object
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/CPCMSolver.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/IEFSolver.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SolverImpl.cpp
  )

# List of headers
list(APPEND headers_list
  CPCMSolver.hpp
  IEFSolver.hpp
  ISolver.hpp
  Solver.hpp
  SolverData.hpp
  SolverImpl.hpp
  )
# Sets install directory for all the headers in the list
foreach(_header IN LISTS headers_list)
  install(FILES ${_header} DESTINATION ${INSTALL_INCLUDEDIR}/${PROJECT_NAME}/solver)
endforeach()
