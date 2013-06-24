cmake_minimum_required(VERSION 2.8)

if( TESTS_NAMES )
	SET( SRC_LIST ${TESTS_NAMES} )
	MESSAGE( STATUS "Building tests: " ${SRC_LIST} )
else()
	execute_process( COMMAND ../updateCMakeSources.sh WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} )
	file(STRINGS CMakeSources.txt SRC_LIST)
	MESSAGE( STATUS "Building all tests" )
endif()

include( ExternalProject )

# Add gtest
ExternalProject_Add(
	googletest
	URL https://googletest.googlecode.com/files/gtest-1.6.0.zip
	UPDATE_COMMAND ""
	INSTALL_COMMAND ""
	)

ExternalProject_Get_Property(googletest source_dir)
include_directories( ../${CMAKE_PROJECT_NAME} ${source_dir}/include)

ExternalProject_Get_Property(googletest binary_dir)
link_directories( ${binary_dir} )

add_executable(test${CMAKE_PROJECT_NAME} ${SRC_LIST})
add_test( ${CMAKE_PROJECT_NAME} test${CMAKE_PROJECT_NAME} )

add_dependencies(test${CMAKE_PROJECT_NAME} googletest)
target_link_libraries(test${CMAKE_PROJECT_NAME} ${EXTERNAL_LIBRARIES} gtest gtest_main pthread )
