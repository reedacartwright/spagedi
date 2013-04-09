# Portable Release Engineering Script

set(PROJECT_VERSION "1.4b-CUSTOM") # Default Version
set(PROJECT_NAME "spagedi")
set(PROJECT_DISTS "SPAGeDi-1*")
set(GITHUB_USER "reedacartwright")
set(GITHUB_PROJECT "${PROJECT_NAME}")
set(GITHUB_URL "http://nodeload.github.com/${GITHUB_USER}/${GITHUB_PROJECT}/legacy.tar.gz/")

# Script Options
# RELENG_TAG
# RELENG_M32
# RELENG_M64
# RELENG_TOOLCHAIN

if(NOT RELENG_TAG)
	set(RELENG_TAG "HEAD")
endif()

string(RANDOM TMP)
set(RELENG_DIR "${CMAKE_CURRENT_BINARY_DIR}/SPAGeDi-temp/${TMP}/")

message(STATUS "Using ${RELENG_DIR} to build packages ...")
file(MAKE_DIRECTORY "${RELENG_DIR}")

set(ARCHIVE "${PROJECT_NAME}.tar.gz")
message(STATUS "Downloading SPAGeDi Archive '${RELENG_TAG}' from Github...")
file(DOWNLOAD "${GITHUB_URL}${RELENG_TAG}" "${RELENG_DIR}${ARCHIVE}" LOG GIT_LOG STATUS GIT_STATUS SHOW_PROGRESS)
list(GET GIT_STATUS 0 GIT_STATUS_0)
list(GET GIT_STATUS 1 GIT_STATUS_1)
if(NOT EXISTS "${RELENG_DIR}${ARCHIVE}" )
	message(ERROR "Download Failed: ${GIT_STATUS_1}")
endif()

# Extract version tag from the download log
if(GIT_LOG MATCHES "Content-Disposition: attachment; filename=${GITHUB_USER}-${GITHUB_PROJECT}-(.*).tar.gz")
	set(PROJECT_VERSION "${CMAKE_MATCH_1}")
	# trim comment id if version matches a tag
	if(PROJECT_VERSION MATCHES "^(.*)-0-.*$")
		set(PROJECT_VERSION "${CMAKE_MATCH_1}")
	endif()
endif()

message(STATUS "Extracting archive ...")
execute_process(COMMAND ${CMAKE_COMMAND} -E tar xzf "${ARCHIVE}"
	WORKING_DIRECTORY "${RELENG_DIR}")
file(GLOB ARCHIVE_DIR "${RELENG_DIR}${GITHUB_USER}-${GITHUB_PROJECT}*")

message(STATUS "Configuring SPAGeDi version ${PROJECT_VERSION} ...")
set(CMAKE_DEFS
	-DCMAKE_BUILD_TYPE=Release
	-DUSE_STATIC_LIBS=on
	-DSPAGEDI_VERSION_GIT=${PROJECT_VERSION}
)
set(CMAKE_ARGS "")
if(WIN32 AND NOT UNIX)
	set(CMAKE_ARGS -G "NMake Makefiles" ${CMAKE_ARGS})
	set(MAKE_BIN nmake)
elseif(APPLE)
	SET(CMAKE_DEFS ${CMAKE_DEFS} 
		"-DCMAKE_OSX_ARCHITECTURES=x86_64\\;i386"
		-DCMAKE_OSX_DEPLOYMENT_TARGET=10.5
		-DCPACK_SYSTEM_NAME=Darwin64-universal
	)
	set(MAKE_BIN make)
else()
	set(MAKE_BIN make)	
endif()

if(RELENG_M32)
	set(CMAKE_DEFS ${CMAKE_DEFS} -DCMAKE_C_FLAGS=-m32)
endif()
if(RELENG_M64)
	set(CMAKE_DEFS ${CMAKE_DEFS} -DCMAKE_C_FLAGS=-m64)
endif()
if(RELENG_TOOLCHAIN)
	get_filename_component(RELENG_TOOLCHAIN ${RELENG_TOOLCHAIN} REALPATH)
	set(CMAKE_DEFS ${CMAKE_DEFS} -DCMAKE_TOOLCHAIN_FILE=${RELENG_TOOLCHAIN})
endif()

execute_process(COMMAND ${CMAKE_COMMAND} ${CMAKE_ARGS} .. ${CMAKE_DEFS}
	WORKING_DIRECTORY "${ARCHIVE_DIR}/build")

message(STATUS "Building packages ...")

execute_process(COMMAND ${MAKE_BIN} ${PROJECT_NAME} WORKING_DIRECTORY "${ARCHIVE_DIR}/build")	
execute_process(COMMAND ${MAKE_BIN} new_package WORKING_DIRECTORY "${ARCHIVE_DIR}/build")	
execute_process(COMMAND ${MAKE_BIN} new_package_source WORKING_DIRECTORY "${ARCHIVE_DIR}/build")	

message(STATUS "Relocating packages ...")
file(GLOB DISTS "${ARCHIVE_DIR}/build/${PROJECT_DISTS}")
foreach(dist ${DISTS})
	execute_process(COMMAND ${CMAKE_COMMAND} -E copy_if_different ${dist} ${CMAKE_CURRENT_BINARY_DIR})
endforeach()

if(NOT NO_CLEAN)
	message(STATUS "Cleaning up ...")
	execute_process(COMMAND ${CMAKE_COMMAND} -E remove_directory "${RELENG_DIR}")
endif()

#get_cmake_property(_variableNames VARIABLES)
#foreach (_variableName ${_variableNames})
#    message(STATUS "${_variableName}=${${_variableName}}")
#endforeach()

