project(SBGS)
cmake_minimum_required(VERSION 2.8)

include_directories( /opt/local/include )

set( CMAKE_CXX_FLAGS "-Wall -Wextra -Werror " )

add_subdirectory( src )
add_subdirectory( tests )

