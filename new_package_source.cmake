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

MESSAGE("Creating Source Package for ${PACKAGE} ${VERSION}...")

foreach(PKG ${CPACK_SOURCE_GENERATOR})
	UNSET(EXT)
	if(PKG STREQUAL "TBZ2")
		SET(EXT "tar.bz2")
	elseif(PKG STREQUAL "TGZ")
		SET(EXT "tar.gz")
	elseif(PKG STREQUAL "ZIP")
		SET(EXT "zip")
	else()
		MESSAGE("ERROR: new_package_source does not understand CPACK_SOURCE_GENERATOR $PKG")
	endif()
	if(DEFINED EXT)
		EXECUTE_PROCESS(COMMAND ${CMAKE_CPACK_COMMAND} -G ${PKG} --config CPackSourceConfig.cmake
			TIMEOUT 3600
			WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

		MESSAGE("Moving ${CPACK_SOURCE_PACKAGE_FILE_NAME}.${EXT} to ${PACKAGE}-${VERSION}.${EXT}")

		EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E rename ${CPACK_SOURCE_PACKAGE_FILE_NAME}.${EXT}
			"${PACKAGE}-${VERSION}.${EXT}")
	endif()
endforeach()


