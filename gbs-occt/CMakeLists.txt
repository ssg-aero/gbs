cmake_minimum_required(VERSION 3.17.0)

project(gbs-occt LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)

option(BUILD_GBS_OCCT_TESTS "Enable the build of gbs-occt C++ tests" OFF)

set(GBS_OCCT_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)

set(GBS_OCCT_HEADERS
    ${GBS_OCCT_INCLUDE_DIR}/gbs/occt/brepbuilders.h
    ${GBS_OCCT_INCLUDE_DIR}/gbs/occt/containers.h
    ${GBS_OCCT_INCLUDE_DIR}/gbs/occt/curvesbuild.h
    ${GBS_OCCT_INCLUDE_DIR}/gbs/occt/discetize.h
    ${GBS_OCCT_INCLUDE_DIR}/gbs/occt/export.h
    ${GBS_OCCT_INCLUDE_DIR}/gbs/occt/funtionsforsolver.h
    ${GBS_OCCT_INCLUDE_DIR}/gbs/occt/geomprim.h
    ${GBS_OCCT_INCLUDE_DIR}/gbs/occt/intersections.h
    ${GBS_OCCT_INCLUDE_DIR}/gbs/occt/measure.h
    ${GBS_OCCT_INCLUDE_DIR}/gbs/occt/surfacesbuild.h
)

set(GBS_OCCT_SRCS
    ${CMAKE_CURRENT_SOURCE_DIR}/src/curvesbuild.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/funtionsforsolver.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/intersections.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/measure.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/surfacesbuild.cpp
)

# occt
find_package(OpenCASCADE REQUIRED)
# gbs
if(NOT TARGET gbs)
    find_package(gbs REQUIRED)
endif ()

add_library(gbs-occt SHARED ${GBS_OCCT_HEADERS} ${GBS_OCCT_SRCS})
add_library(gbs::occt ALIAS gbs-occt)
set_property(TARGET gbs-occt PROPERTY EXPORT_NAME gbs::occt)

target_include_directories(
    gbs-occt
    PUBLIC
        $<BUILD_INTERFACE:${GBS_OCCT_INCLUDE_DIR}>
        $<INSTALL_INTERFACE:include>
        ${OpenCASCADE_INCLUDE_DIR}
)

target_link_libraries(
    gbs-occt
    PUBLIC
        gbs::core
        TKG3d TKG2d TKGeomAlgo TKBRep TKTopAlgo TKOffset TKIGES
)

install(FILES ${GBS_OCCT_HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/gbs/occt)
install(TARGETS gbs-occt)

# C++ tests
if(${BUILD_GBS_OCCT_TESTS})
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
