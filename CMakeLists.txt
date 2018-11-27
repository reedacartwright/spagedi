CMAKE_MINIMUM_REQUIRED(VERSION 2.8 FATAL_ERROR)
SET(CMAKE_VERBOSE_MAKEFILE OFF)
SET(CMAKE_BUILD_TYPE_INIT "RelWithDebInfo")
SET(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/Modules")

PROJECT(SPAGeDi)
SET(SPAGEDI_NAME "SPAGeDi")
SET(SPAGEDI_VERSION "1.5d") # Default Version Number
SET(SPAGEDI_VERSION_MAJOR "1")
SET(SPAGEDI_VERSION_MINOR "5")

SET(SPAGEDI_DIR "spagedi-${SPAGEDI_VERSION_MAJOR}.${SPAGEDI_VERSION_MINOR}")

SET(CMAKE_DATA_DIR "share/${SPAGEDI_DIR}" CACHE STRING
  "Install location for data (relative to prefix).")
SET(CMAKE_DOC_DIR "doc/${SPAGEDI_DIR}" CACHE STRING
  "Install location for documentation (relative to prefix).")
SET(CMAKE_MAN_DIR "man" CACHE STRING
  "Install location for man pages (relative to prefix).")
MARK_AS_ADVANCED(CMAKE_DATA_DIR CMAKE_DOC_DIR CMAKE_MAN_DIR)

ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(example_data)

INSTALL(FILES manual.pdf copying.txt readme.txt DESTINATION ${CMAKE_DOC_DIR})

# uninstall target
if(NOT TARGET uninstall)
    configure_file(
        "${CMAKE_CURRENT_SOURCE_DIR}/cmake_uninstall.cmake.in"
        "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
        IMMEDIATE @ONLY)

    add_custom_target(uninstall
        COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
endif()
