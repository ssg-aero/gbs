cmake_minimum_required(VERSION 3.17.0)

project(gbs-occt LANGUAGES CXX)
set(CMAKE_INCLUDE_CURRENT_DIR ON)
# gbslib
include_directories(${CMAKE_SOURCE_DIR})

file(GLOB OCCT_UTILS_SRC_LIST "*.cpp")
add_library(${PROJECT_NAME} SHARED ${OCCT_UTILS_SRC_LIST})
target_link_libraries(${PROJECT_NAME} 
    ${OpenCASCADE_LIBRARIES}
    TBB::tbb    
)

install(TARGETS ${PROJECT_NAME})
INSTALL_HEADERS(${PROJECT_NAME})
if ( MSVC )
        SET_OUT_DIR(${PROJECT_NAME} ${CMAKE_BINARY_DIR})
endif ( MSVC )
# tests
# add current dir for include in tests
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
if(${GBS_BUILD_TESTS})
    add_subdirectory(tests)
endif(${GBS_BUILD_TESTS})

