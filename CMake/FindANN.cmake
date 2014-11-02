###############################################################################
# Find Flann
#
# This sets the following variables:
# FLANN_FOUND - True if FLANN was found.
# FLANN_INCLUDE_DIRS - Directories containing the FLANN include files.
# FLANN_LIBRARIES - Libraries needed to use FLANN.
# FLANN_DEFINITIONS - Compiler flags for FLANN.
# If FLANN_USE_STATIC is specified and then look for static libraries ONLY else
# look for shared ones

if(ANN_USE_STATIC)
  set(ANN_RELEASE_NAME ANN)
  set(ANN_DEBUG_NAME ANN)
endif(ANN_USE_STATIC)

find_package(PkgConfig)
pkg_check_modules(PC_ANN ANN)
set(ANN_DEFINITIONS ${PC_ANN_CFLAGS_OTHER})

find_path(ANN_INCLUDE_DIR ANN/ANN.h)
		  
find_library(ANN_LIBRARY NAMES ${ANN_RELEASE_NAME})

set(ANN_INCLUDE_DIRS ${ANN_INCLUDE_DIR})
set(ANN_LIBRARIES ${ANN_LIBRARY})

if(ANN_FOUND)
  message(STATUS "ANN found (include: ${ANN_INCLUDE_DIRS}, lib: ${ANN_LIBRARIES})")
endif(ANN_FOUND)