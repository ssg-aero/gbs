cmake_minimum_required(VERSION 3.17.0)

project(render LANGUAGES CXX)
# 
# gbslib
#
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../)
# 
# set files fo lib
# 
file(GLOB GBS_RENDER_SRC_LIST "*.cpp")
message( ${GBS_RENDER_SRC_LIST})
add_library(gbs-${PROJECT_NAME} SHARED ${GBS_RENDER_SRC_LIST})
target_link_libraries(gbs-${PROJECT_NAME} 
    ${VTK_LIBRARIES}
)
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

vtk_module_autoinit(
        TARGETS gbs-${PROJECT_NAME}
        MODULES ${VTK_LIBRARIES}
)

SET_INSTALL_DIR(gbs-${PROJECT_NAME})
