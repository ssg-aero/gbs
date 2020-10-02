cmake_minimum_required(VERSION 3.17.0)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

project(gbslib LANGUAGES CXX)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

if(${GBS_BUILD_TESTS} OR ${USE_OCCT_UTILS})
    add_subdirectory(occt-utils)
endif(${GBS_BUILD_TESTS} OR ${USE_OCCT_UTILS})

if(${USE_RENDER})
    add_subdirectory(render)
endif(${USE_RENDER})

if(${GBS_BUILD_TESTS})
    add_subdirectory(tests)
endif(${GBS_BUILD_TESTS})