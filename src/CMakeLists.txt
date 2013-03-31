
execute_process( COMMAND ../updateCMakeSources.sh WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} )
file(STRINGS CMakeSources.txt SRC_LIST)

add_library(${CMAKE_PROJECT_NAME} ${SRC_LIST})



