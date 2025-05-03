
include_directories(../src)

if (BOOST_SERIALIZATION)
    add_executable(LoadFrictionProblem LoadFrictionProblem.cpp)
    target_link_libraries(LoadFrictionProblem PUBLIC bogus)
endif ()

find_package(fclib CONFIG)
if (fclib_FOUND)
    add_executable(FCLibLoader FCLibLoader.cpp)
    target_link_libraries(FCLibLoader ${PROJECT_NAME} FCLIB::fclib)
endif ()

