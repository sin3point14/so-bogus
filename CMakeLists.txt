project(bogus)
cmake_minimum_required(VERSION 2.6)

OPTION( TESTS "Build tests" ON )
OPTION( APPS "Build apps" ON )
OPTION( LIB "Build lib" ON )

OPTION( OPENMP "Use OpenMP" ON )
OPTION( BOOST_SERIALIZATION "Use boost::serialization" ON )
OPTION( MKL "Use mkl" OFF )

OPTION( WITH_2D "Build lib with 2d support" ON )
OPTION( WITH_3D "Build lib with 3d support" ON )
OPTION( WITH_DYNAMIC "Build lib with dynamic dim support" OFF )

SET( BOGUS_ROOT ${CMAKE_CURRENT_SOURCE_DIR} )
SET( CMAKE_MODULE_PATH ${BOGUS_ROOT}/cmake )


# Eigen
FIND_PACKAGE(Eigen3 REQUIRED)
SET( EXTERNAL_INCLUDES ${EIGEN3_INCLUDE_DIR} )
SET( EXTERNAL_LIBRARIES "" )

# OpenMP
if( OPENMP )
	FIND_PACKAGE(OpenMP)
    if(OPENMP_FOUND)
		set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
		add_definitions(-DEIGEN_DONT_PARALLELIZE)
	else()
    	set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unknown-pragmas " )
    endif()
endif()

# boost
if( BOOST_SERIALIZATION )
                find_package(Boost 1.36.0 COMPONENTS serialization)
                if( Boost_FOUND )
                                add_definitions( -DBOGUS_WITH_BOOST_SERIALIZATION )
                                SET( EXTERNAL_INCLUDES ${EXTERNAL_INCLUDES} ${Boost_INCLUDE_DIR} )
                                SET( EXTERNAL_LIBRARIES ${EXTERNAL_LIBRARIES} boost_serialization )
                endif()
endif()

# MKL
if( MKL )
        find_package( MKL )
        if( MKL_FOUND )
                                add_definitions( -DBOGUS_WITH_MKL )
                                SET( EXTERNAL_INCLUDES ${EXTERNAL_INCLUDES} ${MKL_INCLUDE_DIRS} )
                                SET( EXTERNAL_LIBRARIES ${EXTERNAL_LIBRARIES} ${MKL_LIBRARIES} mkl_intel_lp64 mkl_intel_thread mkl_core iomp5 )
        endif()
endif()

# Build types and processor flags

if (NOT CMAKE_BUILD_TYPE)
  message(STATUS "No build type selected, default to Release")
  set(CMAKE_BUILD_TYPE "Release")
endif()

if ( "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" )
  set (CLANG TRUE)
elseif( "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" )
  set (GCC TRUE)
endif()

if ( GCC OR CLANG )
                SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra " )
                SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Werror " )

                if( GCC AND APPLE )
                                SET( EXTERNAL_LIBRARIES ${EXTERNAL_LIBRARIES} ${EXTERNAL_LIBRARIES} stdc++ )
                endif()

                #Do not use -march=native on gcc < 4.3
                SET(HANDLE_MARCH_NATIVE TRUE)
                if (GCC)
                                execute_process(COMMAND ${CMAKE_C_COMPILER} -dumpversion OUTPUT_VARIABLE GCC_VERSION)
                                if ( NOT GCC_VERSION VERSION_GREATER 4.3 AND NOT GCC_VERSION VERSION_EQUAL 4.3 )
                                                SET(HANDLE_MARCH_NATIVE FALSE)
                                endif()
                endif()
                if( HANDLE_MARCH_NATIVE )
                                SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}  -march=native -DEIGEN_NO_DEBUG ")
                else()
                                SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}  -msse -msse2 -DEIGEN_NO_DEBUG ")
                endif()
endif()

if( MSVC )
        add_definitions( /DBOGUS_DONT_USE_BUILTIN_ATOMICS )
        add_definitions( /wd4305 )
endif()

include_directories( ${EXTERNAL_INCLUDES} )

if( LIB )
        add_subdirectory( src )
endif()
if( TESTS )
                enable_testing()
                add_subdirectory( tests )
endif()
if( APPS )
                add_subdirectory( apps )
endif()
