cmake_minimum_required(VERSION 3.17.0)

project(gbs-mesh LANGUAGES CXX)

set(GBS_MESH_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)

set(GBS_MESH_HEADERS
    ${GBS_MESH_INCLUDE_DIR}/gbs/mesh/mshedge.h
    ${GBS_MESH_INCLUDE_DIR}/gbs/mesh/smoothing.h
    ${GBS_MESH_INCLUDE_DIR}/gbs/mesh/tfi.h
)
add_library(gbs-mesh INTERFACE ${GBS_MESH_HEADERS})
add_library(gbs::mesh ALIAS gbs-mesh)
set_property(TARGET gbs-mesh PROPERTY EXPORT_NAME gbs::mesh)

# gbs
if(NOT TARGET gbs)
    find_package(gbs REQUIRED)
endif ()

target_include_directories(
    ${PROJECT_NAME}
    INTERFACE
        $<BUILD_INTERFACE:${GBS_MESH_INCLUDE_DIR}>
        $<INSTALL_INTERFACE:include>
)
target_link_libraries(
    ${PROJECT_NAME}
    INTERFACE
        gbs::core   
)

install(FILES ${GBS_MESH_HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/gbs/mesh)

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
