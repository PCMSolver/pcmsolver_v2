target_sources(pcm
  PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}/Element.cpp
    ${CMAKE_CURRENT_LIST_DIR}/GePolCavity.cpp
    ${CMAKE_CURRENT_LIST_DIR}/ICavity.cpp
    ${CMAKE_CURRENT_LIST_DIR}/RestartCavity.cpp
  )

# List of headers
list(APPEND headers_list
  Cavity.hpp
  CavityData.hpp
  Element.hpp
  GePolCavity.hpp
  ICavity.hpp
  RestartCavity.hpp
  )
# Sets install directory for all the headers in the list
foreach(_header IN LISTS headers_list)
  install(FILES ${_header} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/cavity)
endforeach()
