cmake_minimum_required(VERSION 3.17.0)

project(gbslib_tests LANGUAGES CXX)
#
# c++ conf
#
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)
set(CMAKE_INCLUDE_CURRENT_DIR ON)
#
# gslib
#
message(${CMAKE_CURRENT_SOURCE_DIR}/../)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../)
#
# CODAN ENV
#
file(TO_CMAKE_PATH $ENV{CONDA_PREFIX} CONDA_PREFIX)
include_directories(${CONDA_PREFIX}/Library/include)
link_directories(${CONDA_PREFIX}/Library/lib)
#
# Locate GTest
#
enable_testing()
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})
#
# eigen
#
include_directories(${CONDA_PREFIX}/Library/include/eigen3)
#
# occt
# 
set(OpenCASCADE_DIR ${CONDA_PREFIX}/Library/lib/cmake/opencascade)
find_package(OpenCASCADE REQUIRED)
link_directories(${OpenCASCADE_LIBRARY_DIR})
include_directories(${OpenCASCADE_INCLUDE_DIR})
add_compile_definitions(_SILENCE_CXX17_ITERATOR_BASE_CLASS_DEPRECATION_WARNING) # occt interators seems to be old
#
# occt-utils
#
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../occt-utils)

file(GLOB SRC_LIST "tests_*")
message(${SRC_LIST})
add_executable(gbslib_tests ${SRC_LIST})

target_link_libraries(gbslib_tests 
                            GTest::GTest 
                            GTest::Main 
                            ${OpenCASCADE_LIBRARIES} 
                            occt-utils
                            sundials_kinsol
                            nlopt
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