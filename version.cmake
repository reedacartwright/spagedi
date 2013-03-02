#get_cmake_property(_variableNames VARIABLES)
#foreach (_variableName ${_variableNames})
#    message(STATUS "${_variableName}=${${_variableName}}")
#endforeach()

SET(VERFILE "src/spagedi_version.h")

# If VERFILE exists in src, copy it to build
IF(EXISTS "${SRC}/${VERFILE}")
	EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E copy
		"${SRC}/${VERFILE}"
	    "${DST}/${VERFILE}")
ENDIF()

# Update version if we can
IF(GIT_EXECUTABLE AND EXISTS "${SRC}/.git") 
	EXECUTE_PROCESS(
	     COMMAND ${GIT_EXECUTABLE} describe --tags
	     OUTPUT_VARIABLE VERSION
	     OUTPUT_STRIP_TRAILING_WHITESPACE
	)
	CONFIGURE_FILE(${SRC}/${VERFILE}.in
                   ${DST}/${VERFILE}
                   @ONLY)
ENDIF()

# Use Default version if it doesn't exist
IF(NOT EXISTS "${DST}/${VERFILE}")
	SET(VERSION ${VER})
	CONFIGURE_FILE(${SRC}/${VERFILE}.in
                   ${DST}/${VERFILE}
                   @ONLY)	
ENDIF()

