# Choose an appropriate compiler prefix
SET(COMPILER_PREFIX i486-mingw32) # Arch Linux

# the name of the target operating system
SET(CMAKE_SYSTEM_NAME Windows)
SET(CMAKE_SYSTEM_PROCESSOR x86)
 
# which compilers to use for C and C++
find_program(CMAKE_RC_COMPILER  NAMES ${COMPILER_PREFIX}-windres)
find_program(CMAKE_C_COMPILER   NAMES ${COMPILER_PREFIX}-gcc)
find_program(CMAKE_CXX_COMPILER NAMES ${COMPILER_PREFIX}-g++)

# here is the target environment located
#SET(USER_ROOT_PATH)
SET(CMAKE_FIND_ROOT_PATH /usr/${CROSS_PREFIX} /usr/local/${CROSS_PREFIX})
 
# adjust the default behaviour of the FIND_XXX() commands:
# search headers and libraries in the target environment, search 
# programs in the host environment
set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)

