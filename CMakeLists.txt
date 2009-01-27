CMAKE_MINIMUM_REQUIRED(VERSION 2.6 FATAL_ERROR)
PROJECT(SPAGeDi)

SET(SPAGEDI_VERSION "1.2i")

INCLUDE(CheckIncludeFile)
INCLUDE(CheckFunctionExists)
INCLUDE(CheckLibraryExists)

SET(CMAKE_VERBOSE_MAKEFILE ON)

SET(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/Modules")

SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "SPAGeDi")
SET(CPACK_PACKAGE_VENDOR "Olivier Hardy & Xavier Vekemans")
SET(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/readme.txt")
SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/copying.txt")
SET(CPACK_PACKAGE_VERSION "${SPAGEDI_VERSION}")

#Bug in cpack
SET(CPACK_SYSTEM_NAME ${CMAKE_SYSTEM_NAME})
IF(${CPACK_SYSTEM_NAME} MATCHES Windows)
  IF(CMAKE_CL_64)
    SET(CPACK_SYSTEM_NAME win64)
  ELSE(CMAKE_CL_64)
    SET(CPACK_SYSTEM_NAME win32)
  ENDIF(CMAKE_CL_64)
ENDIF(${CPACK_SYSTEM_NAME} MATCHES Windows)

SET(CPACK_PACKAGE_FILE_NAME "${CMAKE_PROJECT_NAME}-${CPACK_PACKAGE_VERSION}-${CPACK_SYSTEM_NAME}")
SET(CPACK_SOURCE_PACKAGE_FILE_NAME "${CMAKE_PROJECT_NAME}-${CPACK_PACKAGE_VERSION}-Source")

SET(CPACK_PACKAGE_INSTALL_DIRECTORY "SPAGeDi")

SET(CPACK_SOURCE_IGNORE_FILES
  "/CVS/"
  "/\\\\.svn/"
  "\\\\.swp$"
  "\\\\.#"
  "/#"
  ".*~$"
  "SPAGeDi[-]"
  "spagedi\\\\.exe$"
  "spagedi$"
  "/CMakeFiles/"
  "CMakeCache\\\\.txt"
  "CPack.*Config\\\\.cmake"
  "cmake_install\\\\.cmake"
  "_CPACK_PACKAGES"
)

IF(WIN32 AND NOT UNIX)
  # There is a bug in NSI that does not handle full unix paths properly. Make
  # sure there is at least one set of four (4) backlasshes.
  #SET(CPACK_PACKAGE_ICON "${CMake_SOURCE_DIR}/Utilities/Release\\\\InstallIcon.bmp")
  SET(CPACK_NSIS_INSTALLED_ICON_NAME "bin\\\\spagedi.exe")
  SET(CPACK_NSIS_DISPLAY_NAME "${CPACK_PACKAGE_INSTALL_DIRECTORY} SPAGeDi")
  SET(CPACK_NSIS_HELP_LINK "http:\\\\\\\\www.ulb.ac.be\\\\sciences\\\\ecoevol\\\\spagedi.html")
  SET(CPACK_NSIS_URL_INFO_ABOUT "http:\\\\\\\\www.ulb.ac.be\\\\sciences\\\\ecoevol\\\\spagedi.html")
  SET(CPACK_NSIS_CONTACT "ohardy@ulb.ac.be")
  SET(CPACK_NSIS_MODIFY_PATH ON)
ELSE(WIN32 AND NOT UNIX)
  SET(CPACK_GENERATOR "TGZ")
  SET(CPACK_SOURCE_GENERATOR "TGZ") 
  SET(CPACK_STRIP_FILES "bin/spagedi")
ENDIF(WIN32 AND NOT UNIX)

SET(CPACK_PACKAGE_EXECUTABLES "spagedi" "SPAGeDi")
INCLUDE(CPack)

INSTALL(FILES manual.pdf copying.txt readme.txt DESTINATION share)
INSTALL(FILES examples/data1.txt examples/data2.txt examples/data3.txt DESTINATION share/examples)


ADD_SUBDIRECTORY(src)

INCLUDE(InstallRequiredSystemLibraries)
 
