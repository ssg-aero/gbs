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

# gbs
if(NOT TARGET gbs)
    find_package(gbs REQUIRED)
endif ()

target_include_directories(
    gbs-mesh
    INTERFACE
        $<BUILD_INTERFACE:${GBS_MESH_INCLUDE_DIR}>
        $<INSTALL_INTERFACE:include>
        gbs::core
)

install(FILES ${GBS_MESH_HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/gbs/mesh)
