cmake_minimum_required(VERSION 3.17.0)

project(gbs_bindings LANGUAGES CXX)

# pybind11
find_package(pybind11 REQUIRED)
# gbs
if(NOT TARGET gbs)
    find_package(gbs REQUIRED)
endif ()

set(PYGBS_SRCS
    ${CMAKE_CURRENT_SOURCE_DIR}/src/gbsbindrender.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/gbsbindmesh.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/gbsbindinterp.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/gbsbind.cpp
)

pybind11_add_module(
    ${PROJECT_NAME}
    ${PYGBS_SRCS} 
)
target_link_libraries(
    ${PROJECT_NAME}
    PRIVATE 
        gbs::core
        gbs::render
        gbs::mesh
)

set_target_properties(
    ${PROJECT_NAME} PROPERTIES
    OUTPUT_NAME "gbs")

install(TARGETS ${PROJECT_NAME} 
        LIBRARY DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/pygbs/)
