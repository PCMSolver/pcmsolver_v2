if(NOT DEFINED PYMOD_INSTALL_LIBDIR)
  message(STATUS "Setting (unspecified) option PYMOD_INSTALL_LIBDIR: python")
  set(PYMOD_INSTALL_LIBDIR "python" CACHE STRING "Location within INSTALL_LIBDIR to which Python modules are installed" FORCE)
else()
  message(STATUS "Setting option PYMOD_INSTALL_LIBDIR: ${PYMOD_INSTALL_LIBDIR}")
  set(PYMOD_INSTALL_LIBDIR "${PYMOD_INSTALL_LIBDIR}" CACHE STRING "Location within INSTALL_LIBDIR to which Python modules are installed" FORCE)
endif()
file(TO_NATIVE_PATH "${INSTALL_LIBDIR}/${PYMOD_INSTALL_LIBDIR}/pcmsolver" PYMOD_INSTALL_FULLDIR)

configure_file(go_pcm.py.in ${PROJECT_BINARY_DIR}/${INSTALL_BINDIR}/go_pcm.py @ONLY)
file(
  COPY
    plot_cavity.py
  DESTINATION
    ${PROJECT_BINARY_DIR}/${INSTALL_BINDIR}
  )
install(
  FILES
    ${PROJECT_BINARY_DIR}/${INSTALL_BINDIR}/go_pcm.py
    ${PROJECT_BINARY_DIR}/${INSTALL_BINDIR}/plot_cavity.py
  DESTINATION
    ${INSTALL_BINDIR}
  )

configure_file(init.py ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}/__init__.py COPYONLY)
file(
  COPY
    getkw.py
    pcmdata.py
    pcmparser.py
    pyparsing.py
  DESTINATION
    ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}
  )
install(
  FILES
    ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}/__init__.py
    ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}/getkw.py
    ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}/pcmdata.py
    ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}/pcmparser.py
    ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}/pyparsing.py
  DESTINATION
    ${PYMOD_INSTALL_FULLDIR}
  )
