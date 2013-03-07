# Command for creating a package source with defined name

INCLUDE(CPackSourceConfig.cmake)

#get_cmake_property(_variableNames VARIABLES)
#foreach (_variableName ${_variableNames})
#    message(STATUS "${_variableName}=${${_variableName}}")
#endforeach()

FILE(STRINGS ${VERFILE} VERSION LIMIT_COUNT 1)

if("${VERSION}" MATCHES "\"(.*) +(.*)\"")
	SET(VERSION "${CMAKE_MATCH_2}")
	SET(PACKAGE "${CMAKE_MATCH_1}")
else()
	SET(VERSION "${CPACK_PACKAGE_VERSION})
	SET(PACKAGE "${CPACK_PACKAGE_NAME})
endif()

#MESSAGE("Creating Source Package for ${PACKAGE} ${VERSION}...")

EXECUTE_PROCESS(COMMAND ${CMAKE_CPACK_COMMAND}
	-D CPACK_PACKAGE_FILE_NAME=${PACKAGE}-${VERSION}
	-D VERFILE=${VERFILE}
	--config CPackSourceConfig.cmake
	-P "${PACKAGE}"
	-R "${VERSION}"
	TIMEOUT 3600
	WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

