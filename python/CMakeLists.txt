cmake_minimum_required(VERSION 3.17.0)

project(gbs LANGUAGES CXX)
#
# gslib
#
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../)

find_package(pybind11 REQUIRED)
pybind11_add_module(gbs gbsbind.cpp)
target_link_libraries(gbs PRIVATE nlopt)

if ( MSVC )
    SET_OUT_DIR(${PROJECT_NAME} ${CMAKE_BINARY_DIR})
endif ( MSVC )