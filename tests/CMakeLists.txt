
execute_process( COMMAND ../updateCMakeSources.sh WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} )
file(STRINGS CMakeSources.txt SRC_LIST)

include_directories( ../src )

add_executable(test${CMAKE_PROJECT_NAME} ${SRC_LIST})
target_link_libraries(test${CMAKE_PROJECT_NAME} lib${CMAKE_PROJECT_NAME})

