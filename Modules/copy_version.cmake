# This shouldn't be run if we are using the built in packaging targets
#IF(NOT NEW_PACKAGE)
#	RETURN()
#ENDIF()

#get_cmake_property(_variableNames VARIABLES)
#foreach (_variableName ${_variableNames})
#    message("    ${_variableName}=${${_variableName}}")
#endforeach()

SET(VERFILE_BASE "src/spagedi_version.h")

# Test if we are doing a source install
IF(NOT CPACK_INSTALL_CMAKE_PROJECTS)
	IF(EXISTS "${CPACK_BINARY_DIR}/${VERFILE_BASE}")
		# Copy Version File
		EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E copy
			"${CPACK_BINARY_DIR}/${VERFILE_BASE}"
	        "${CMAKE_INSTALL_PREFIX}/${VERFILE_BASE}.pkg")
	ELSE()
		# Constuct a pkg version file
		SET(PACKAGE_NAME ${CPACK_PACKAGE_NAME})
		SET(PACKAGE_VERSION ${CPACK_PACKAGE_VERSION})
		CONFIGURE_FILE("${CPACK_SOURCE_DIR}/${VERFILE_BASE}.in"
			"${CMAKE_INSTALL_PREFIX}${VERFILE_BASE}.pkg"
			@ONLY)
	ENDIF()
	
	# Construct Build Directory
	FILE(MAKE_DIRECTORY "${CMAKE_INSTALL_PREFIX}/build")
	FILE(WRITE "${CMAKE_INSTALL_PREFIX}/build/build.txt"
		"Read ../readme.txt for building instructions.\n")
ENDIF()

