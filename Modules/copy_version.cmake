# Test if we are doing a source install
IF(NOT CPACK_INSTALL_CMAKE_PROJECTS)
	IF(EXISTS "${VERFILE}")
		# Copy version.h
		EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E copy
			"${VERFILE}"
	        "${CMAKE_CURRENT_BINARY_DIR}/src/spagedi_version.h")
	ENDIF()
	
	# Construct Build Directory
	FILE(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/build")
	FILE(WRITE "${CMAKE_CURRENT_BINARY_DIR}/build/build.txt"
		"Read ../readme.txt for building instructions.\n")
ENDIF()

