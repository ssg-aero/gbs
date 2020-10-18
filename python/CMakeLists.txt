cmake_minimum_required(VERSION 3.17.0)

project(pygbs LANGUAGES CXX)
#
# gslib
#
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../)

find_package(pybind11 REQUIRED)
pybind11_add_module(${PROJECT_NAME} gbsbind.cpp)
target_link_libraries(${PROJECT_NAME} PRIVATE nlopt)

if ( MSVC )
    SET_OUT_DIR(${PROJECT_NAME} ${CMAKE_BINARY_DIR})
endif ( MSVC )