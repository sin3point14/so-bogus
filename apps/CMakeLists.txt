
include_directories( ../src )

if( BOOST_SERIALIZATION )
                add_executable (LoadFrictionProblem LoadFrictionProblem.cpp)
                target_link_libraries (LoadFrictionProblem PUBLIC bogus)
endif()

find_package(FCLib)
if( FCLIB_FOUND )
        include_directories( ${FCLib_INCLUDE_DIRS} )
        add_executable (FCLibLoader FCLibLoader.cpp)
        target_link_libraries (FCLibLoader ${PROJECT_NAME} ${FCLib_LIBRARIES} )
endif()

