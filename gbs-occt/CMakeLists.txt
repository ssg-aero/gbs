cmake_minimum_required(VERSION 3.17.0)

project(occt LANGUAGES CXX)
# 
# gbslib
#
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../)
# 
# set files fo lib
# 
file(GLOB OCCT_UTILS_SRC_LIST "*.cpp")
add_library(gbs-${PROJECT_NAME} SHARED ${OCCT_UTILS_SRC_LIST})
target_link_libraries(gbs-${PROJECT_NAME} ${OpenCASCADE_LIBRARIES})

if ( MSVC )
    SET_OUT_DIR(gbs-${PROJECT_NAME} ${CMAKE_BINARY_DIR})
endif ( MSVC )
#
# tests
#
# add current dir for include in tests
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
if(${GBS_BUILD_TESTS})
    add_subdirectory(tests)
endif(${GBS_BUILD_TESTS})

SET_INSTALL_DIR(gbs-${PROJECT_NAME})