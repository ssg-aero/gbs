cmake_minimum_required(VERSION 3.17.0)

project(pygbs LANGUAGES CXX)
#
# gslib
#
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../)

find_package(pybind11 REQUIRED)
pybind11_add_module(gbs gbsbind.cpp)
target_link_libraries(gbs PRIVATE nlopt)

#
# set output dirs (vs is a pain)
#
if ( MSVC )
    foreach( OUTPUTCONFIG ${CMAKE_CONFIGURATION_TYPES} )
        string( TOUPPER ${OUTPUTCONFIG} OUTPUTCONFIG )
        set_property(TARGET gbs PROPERTY RUNTIME_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${CMAKE_BINARY_DIR} )
        set_property(TARGET gbs PROPERTY LIBRARY_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${CMAKE_BINARY_DIR} )
        set_property(TARGET gbs PROPERTY ARCHIVE_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${CMAKE_BINARY_DIR} )
    endforeach( OUTPUTCONFIG CMAKE_CONFIGURATION_TYPES )
endif ( MSVC )