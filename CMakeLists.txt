project(bogus)
cmake_minimum_required(VERSION 2.8)

SET( BOGUS_ROOT ${CMAKE_CURRENT_SOURCE_DIR} )
SET( CMAKE_MODULE_PATH ${BOGUS_ROOT}/cmake )
# Eigen
FIND_PACKAGE(Eigen3 REQUIRED)
SET( EXTERNAL_INCLUDES ${EIGEN3_INCLUDE_DIR} )
add_definitions(-DEIGEN_DONT_PARALLELIZE)
SET( EXTERNAL_LIBRARIES "" )

# OpenMP / pthreads
FIND_PACKAGE(OpenMP)
if(OPENMP_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
  SET( EXTERNAL_LIBRARIES ${EXTERNAL_LIBRARIES} ${OpenMP_LIBRARIES} )
else()
set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unknown-pragmas " )
endif()

# Proc flags
set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Werror -g " )

include_directories( ${EXTERNAL_INCLUDES} )

add_subdirectory( src )
add_subdirectory( tests )

