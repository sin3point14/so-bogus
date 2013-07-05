# - Try to find MKL lib
#
#
# Once done this will define
#
#  MKL_FOUND - system has MKL lib
#  MKL_INCLUDE_DIRS - the MKL include directory
#  MKL_LIBRARIES - link this to use MKL

find_library(MKL_LIBRARIES NAMES  mkl_core 
  PATHS 
  ${CMAKE_INSTALL_PREFIX}/lib/intel/em64t
  $ENV{MKL_ROOT}/lib/em64t
  )
find_path(MKL_INCLUDE_DIRS NAMES mkl.h
  PATHS
  ${CMAKE_INSTALL_PREFIX}/include
  $ENV{MKL_ROOT}/include
  )
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MKL "MKL library not found ; set environment variable MKL_ROOT " MKL_LIBRARIES MKL_INCLUDE_DIRS)
#make_library_set( ${MKL_LIBRARIES} )
mark_as_advanced( MKL_LIBRARIES )
mark_as_advanced( MKL_INCLUDE_DIRS )


