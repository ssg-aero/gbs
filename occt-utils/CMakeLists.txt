cmake_minimum_required(VERSION 3.17.0)

project(occt-utils LANGUAGES CXX)
# c++ conf
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)
set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(HOME_DIR $ENV{USERPROFILE})
# 
# gbslib
#
include_directories(${HOME_DIR}/workspace/gbslib)
# CODAN ENV
file(TO_CMAKE_PATH $ENV{CONDA_PREFIX} CONDA_PREFIX)
#
# occt
# 
set(OpenCASCADE_DIR ${CONDA_PREFIX}/Library/lib/cmake/opencascade)
find_package(OpenCASCADE REQUIRED)
link_directories(${OpenCASCADE_LIBRARY_DIR})
message(${OpenCASCADE_LIBRARY_DIR})
include_directories(${OpenCASCADE_INCLUDE_DIR})
add_compile_definitions(_SILENCE_CXX17_ITERATOR_BASE_CLASS_DEPRECATION_WARNING) # occt iterators seems to be old
# 
# set files fo lib
# 
file(GLOB OCCT_UTILS_SRC_LIST "occt-utils/*.cpp")
add_library(${PROJECT_NAME} SHARED ${OCCT_UTILS_SRC_LIST})
target_link_libraries(${PROJECT_NAME} ${OpenCASCADE_LIBRARIES})
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
add_subdirectory(tests)
