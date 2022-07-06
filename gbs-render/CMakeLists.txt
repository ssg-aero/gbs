cmake_minimum_required(VERSION 3.17.0)

project(gbs-render LANGUAGES CXX)

option(BUILD_GBS_RENDER_TESTS "Enable the build of gbs-render C++ tests" OFF)

set(GBS_RENDER_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)

set(GBS_RENDER_HEADERS
    ${GBS_RENDER_INCLUDE_DIR}/gbs/render/vtkcurvesrender.h
    ${GBS_RENDER_INCLUDE_DIR}/gbs/render/vtkfunctionrender.h
    ${GBS_RENDER_INCLUDE_DIR}/gbs/render/vtkgridrender.h
)

set(GBS_RENDER_SRCS
    ${CMAKE_CURRENT_SOURCE_DIR}/src/vtkcurvesrender.cpp
)

# gbs-core
if(NOT TARGET gbs)
    find_package(gbs REQUIRED)
endif ()
# VTK
find_package(VTK COMPONENTS
    CommonCore
    CommonColor
    CommonDataModel
    FiltersSources
    FiltersGeneral
    InteractionStyle
    InteractionWidgets
    RenderingCore
    RenderingFreeType
    RenderingOpenGL2
    RenderingAnnotation
    ChartsCore
    RenderingContext2D
    RenderingContextOpenGL2
    ViewsContext2D
    RenderingGL2PSOpenGL2
    IOXML 
    Python
    WrappingPythonCore
)

add_library(gbs-render SHARED ${GBS_RENDER_HEADERS} ${GBS_RENDER_SRCS})
add_library(gbs::render ALIAS gbs-render)
set_property(TARGET gbs-render PROPERTY EXPORT_NAME gbs::render)

target_include_directories(
    gbs-render
    PUBLIC
        $<BUILD_INTERFACE:${GBS_RENDER_INCLUDE_DIR}>
        $<INSTALL_INTERFACE:include>
)

target_link_libraries(
    gbs-render
    PUBLIC
        gbs::core
        ${VTK_LIBRARIES} 
)

install(FILES ${GBS_RENDER_HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/gbs/render)
install(TARGETS gbs-render)

# C++ tests
if(${BUILD_GBS_RENDER_TESTS})
    add_subdirectory(tests)
endif()

# CMake config
# Add helpers for cmake configuration
include(CMakePackageConfigHelpers)

# Definition
set(ConfigFileInstallDir ${CMAKECONFIG_INSTALL_DIR}/${PROJECT_NAME})
configure_package_config_file(${PROJECT_NAME}Config.cmake.in
                            "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
                            INSTALL_DESTINATION ${ConfigFileInstallDir})

# Makes the project importable from the build directory
install(TARGETS ${PROJECT_NAME} EXPORT ${PROJECT_NAME}-targets)
export(EXPORT ${PROJECT_NAME}-targets
       FILE "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Targets.cmake")

# Installation
install(EXPORT ${PROJECT_NAME}-targets 
        FILE ${PROJECT_NAME}Targets.cmake
        DESTINATION ${ConfigFileInstallDir})

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
        DESTINATION ${ConfigFileInstallDir})
