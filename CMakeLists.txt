cmake_minimum_required(VERSION 3.17.0)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
project(gbslib LANGUAGES CXX)

# add current dir for include in tests
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
add_subdirectory(tests)