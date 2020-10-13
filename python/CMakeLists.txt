cmake_minimum_required(VERSION 3.17.0)

project(pygbs LANGUAGES CXX)
#
# c++ conf
#
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)
set(CMAKE_INCLUDE_CURRENT_DIR ON)
#
# CODAN ENV
#
file(TO_CMAKE_PATH $ENV{CONDA_PREFIX} CONDA_PREFIX)
include_directories(${CONDA_PREFIX}/Library/include)
include_directories(${CONDA_PREFIX}/include)
link_directories(${CONDA_PREFIX}/Library/lib)
link_directories(${CONDA_PREFIX}/Library/bin)
#
# eigen
#
include_directories(${CONDA_PREFIX}/Library/include/eigen3)
#
# gslib
#
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../)

# file(GLOB SRC_LIST "*.cpp")
# message(${SRC_LIST})
# add_library(${PROJECT_NAME} ${SRC_LIST})

# target_link_libraries(${PROJECT_NAME} 

# add_subdirectory(pybind11)
find_package(pybind11 REQUIRED)
pybind11_add_module(gbs gbsbind.cpp)
target_link_libraries(gbs PRIVATE nlopt)
# pybind11_add_module(example example.cpp)

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