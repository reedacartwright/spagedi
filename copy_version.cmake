# Test if we are doing a source install
# And install version.h into tarball
IF(NOT CPACK_INSTALL_CMAKE_PROJECTS)
IF(EXISTS "${VERFILE}")
	EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E copy
		"${VERFILE}"
        "${CMAKE_CURRENT_BINARY_DIR}/src/version.h")
ENDIF()
ENDIF()

