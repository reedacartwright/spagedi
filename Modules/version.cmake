#get_cmake_property(_variableNames VARIABLES)
#foreach (_variableName ${_variableNames})
#    message(STATUS "${_variableName}=${${_variableName}}")
#endforeach()

SET(VERFILE_BASE "src/spagedi_version.h")

# Update version if we can
IF(GIT_EXECUTABLE AND EXISTS "${SRC}/.git") 
	EXECUTE_PROCESS(
	     COMMAND ${GIT_EXECUTABLE} describe --tags --dirty
	     RESULT_VARIABLE GIT_RESULT
	     OUTPUT_VARIABLE PKG_VERSION
	     OUTPUT_STRIP_TRAILING_WHITESPACE
	     ERROR_QUIET
	)
	IF(GIT_RESULT EQUAL 0)
		SET(PACKAGE_VERSION "${PKG_VERSION}")
	ENDIF()
ENDIF()

# Prefer .pkg which
IF(EXISTS ${SRC}/${VERFILE_BASE}.pkg)
	SET(INFILE ${SRC}/${VERFILE_BASE}.pkg)
ELSE()
	SET(INFILE ${SRC}/${VERFILE_BASE}.in)
ENDIF()

CONFIGURE_FILE(${INFILE} ${DST}/${VERFILE_BASE} @ONLY)

