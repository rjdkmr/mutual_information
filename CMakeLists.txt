cmake_minimum_required(VERSION 2.8)
project(mi CXX)

set(mi_MAJOR_VERSION 1)
set(mi_MINOR_VERSION 0)
set(mi_PATCH_VERSION 0)
set(mi_VERSION
  ${mi_MAJOR_VERSION}.${mi_MINOR_VERSION}.${mi_PATCH_VERSION})

# Offer the user the choice of overriding the installation directories
set(INSTALL_LIB_DIR lib CACHE PATH "Installation directory for libraries")
set(INSTALL_BIN_DIR bin CACHE PATH "Installation directory for executables")
set(INSTALL_INCLUDE_DIR include CACHE PATH "Installation directory for header files")
if(WIN32 AND NOT CYGWIN)
	set(DEF_INSTALL_CMAKE_DIR CMake)
else()
set(DEF_INSTALL_CMAKE_DIR lib/CMake/mi)
endif()
set(INSTALL_CMAKE_DIR ${DEF_INSTALL_CMAKE_DIR} CACHE PATH "Installation directory for CMake files")

# Make relative paths absolute (needed later on)
foreach(p LIB BIN INCLUDE CMAKE)
	set(var INSTALL_${p}_DIR)
	if(NOT IS_ABSOLUTE "${${var}}")
		set(${var} "${CMAKE_INSTALL_PREFIX}/${${var}}")
	endif()
endforeach()

# set up include-directories
include_directories(
	"${mi_SOURCE_DIR}"   # to find src/*.h
	"${mi_BINARY_DIR}")  # to find src/config.h
	 
# Add sub-directories
add_subdirectory(src)


# ============== Eclipse CDT Error Handling======
IF(CMAKE_COMPILER_IS_GNUCC)
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fmessage-length=0")
ENDIF(CMAKE_COMPILER_IS_GNUCC)
IF(CMAKE_COMPILER_IS_GNUCXX)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fmessage-length=0")
ENDIF(CMAKE_COMPILER_IS_GNUCXX)
set(CMAKE_VERBOSE_MAKEFILE ON)

# ===============================