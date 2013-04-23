cmake_minimum_required(VERSION 2.8)

execute_process( COMMAND ../updateCMakeSources.sh WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} )
file(STRINGS CMakeSources.txt SRC_LIST)

include_directories( ../src ${source_dir}/include)

FOREACH( test_file ${SRC_LIST} )
  STRING(REGEX REPLACE "([a-zA-Z0-9]+).cpp" "\\1" test_name ${test_file} )

  add_executable (${test_name} ${test_file})
  target_link_libraries (${test_name} ${CMAKE_PROJECT_NAME} )

ENDFOREACH()
