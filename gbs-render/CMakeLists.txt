cmake_minimum_required(VERSION 3.17.0)

project(gbs-render LANGUAGES CXX)
# 
# gbslib
#
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../)
# 
# set files fo lib
# 
file(GLOB GBS_RENDER_SRC_LIST "*.cpp")
message( ${GBS_RENDER_SRC_LIST})
add_library(${PROJECT_NAME} SHARED ${GBS_RENDER_SRC_LIST})
target_link_libraries(${PROJECT_NAME} 
    ${VTK_LIBRARIES}
)
#
# set output dirs (vs is a pain)
#
if ( MSVC )
    foreach( OUTPUTCONFIG ${CMAKE_CONFIGURATION_TYPES} )
        string( TOUPPER ${OUTPUTCONFIG} OUTPUTCONFIG )
        set_property(TARGET ${PROJECT_NAME} PROPERTY RUNTIME_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${CMAKE_BINARY_DIR} )
        set_property(TARGET ${PROJECT_NAME} PROPERTY LIBRARY_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${CMAKE_BINARY_DIR} )
        set_property(TARGET ${PROJECT_NAME} PROPERTY ARCHIVE_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${CMAKE_BINARY_DIR} )
    endforeach( OUTPUTCONFIG CMAKE_CONFIGURATION_TYPES )
endif ( MSVC )
#
# tests
#
# add current dir for include in tests

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
if(${GBS_BUILD_TESTS})
    add_subdirectory(tests)
endif(${GBS_BUILD_TESTS})

vtk_module_autoinit(
        TARGETS ${PROJECT_NAME}
        MODULES ${VTK_LIBRARIES}
)

SET_INSTALL_DIR(${PROJECT_NAME})
