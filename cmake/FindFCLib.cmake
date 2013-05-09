# - Try to find FCLib lib
#
#
# Once done this will define
#
#  FCLIB_FOUND - system has FCLib lib
#  FCLib_INCLUDE_DIRS - the FCLib include directory
#  FCLib_LIBRARIES - link this to use FCLib

find_library(FCLib_LIBRARIES NAMES fclib
  PATHS 
  ${CMAKE_INSTALL_PREFIX}/lib
  $ENV{FCLIB_ROOT}/lib
  )
find_path(FCLib_INCLUDE_DIRS NAMES fclib.h
  PATHS
  ${CMAKE_INSTALL_PREFIX}/include
  $ENV{FCLIB_ROOT}/include
  )
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FCLib "FCLib library not found ; set environment variable FCLIB_ROOT " FCLib_LIBRARIES FCLib_INCLUDE_DIRS)
#make_library_set( ${FCLib_LIBRARIES} )
mark_as_advanced( FCLib_LIBRARIES )
mark_as_advanced( FCLib_INCLUDE_DIRS )


