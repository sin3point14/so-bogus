execute_process( COMMAND ../updateCMakeSources.sh WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} )
file(STRINGS CMakeSources.txt SRC_LIST)

if( WITH_2D )
        add_definitions( -DBOGUS_INSTANTIATE_2D_SOC )
endif()
if( WITH_3D )
        add_definitions( -DBOGUS_INSTANTIATE_3D_SOC )
endif()
if( WITH_DYNAMIC )
        add_definitions( -DBOGUS_INSTANTIATE_DYNAMIC_SOC )
endif()

#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -S --verbose-asm " )

add_library(bogus STATIC ${SRC_LIST})
target_link_libraries(bogus ${EXTERNAL_LIBRARIES} )

FOREACH(item ${SRC_LIST})
  STRING(REGEX MATCH "^bogus/Interfaces/.+\\.hpp$" item ${item})
  IF(item)
	LIST(APPEND INTERFACES ${item})
  ENDIF(item)
ENDFOREACH(item ${SRC_LIST})

INSTALL( DIRECTORY bogus DESTINATION include )
INSTALL( FILES ${INTERFACES} DESTINATION include/bogus/Interfaces )
INSTALL( TARGETS bogus
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION lib
  ARCHIVE DESTINATION lib )
