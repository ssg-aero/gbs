cmake_minimum_required(VERSION 3.17.0)

project(gbs_render LANGUAGES CXX)
# c++ conf
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)
set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(HOME_DIR $ENV{USERPROFILE})
# 
# gbslib
#
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../)
#
# current tested lib
# 
link_directories(${CMAKE_BINARY_DIR})
# CODAN ENV
file(TO_CMAKE_PATH $ENV{CONDA_PREFIX} CONDA_PREFIX)
link_directories(${CONDA_PREFIX}/Library/lib)
# 
# VTK
# 
set(VTK_DIR  ${CONDA_PREFIX}/Library/lib/cmake/vtk-9.0)
include_directories(${CONDA_PREFIX}/Library/include/vtk-9.0)
find_package(VTK COMPONENTS 
CommonCore
CommonColor
CommonDataModel
FiltersSources
InteractionStyle
RenderingCore
RenderingFreeType
RenderingOpenGL2
)
# 
# set files fo lib
# 
file(GLOB GBS_RENDER_SRC_LIST "*.cpp")
message( ${GBS_RENDER_SRC_LIST})
add_library(${PROJECT_NAME} SHARED ${GBS_RENDER_SRC_LIST})
target_link_libraries(${PROJECT_NAME} 
    ${VTK_LIBRARIES}
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
#
# tests
#
# add current dir for include in tests
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
add_subdirectory(tests)

vtk_module_autoinit(
        TARGETS ${PROJECT_NAME}
        MODULES ${VTK_LIBRARIES}
)