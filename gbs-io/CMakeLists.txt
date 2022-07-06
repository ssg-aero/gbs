cmake_minimum_required(VERSION 3.17.0)

project(gbs-io LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)

set(GBS_IO_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)

set(GBS_IO_HEADERS
    ${GBS_IO_INCLUDE_DIR}/gbs/io/fromjson.h
    ${GBS_IO_INCLUDE_DIR}/gbs/io/fromtext.h
    ${GBS_IO_INCLUDE_DIR}/gbs/io/iges.h
    ${GBS_IO_INCLUDE_DIR}/gbs/io/print.h
    ${GBS_IO_INCLUDE_DIR}/gbs/io/tojson.h
)
set(GBS_IO_3RD_PARTIES 
    ${GBS_IO_INCLUDE_DIR}/gbs/third_parties/magic_enum.hpp
)

add_library(gbs-io INTERFACE ${GBS_IO_HEADERS} ${GBS_IO_3RD_PARTIES})
add_library(gbs::io ALIAS gbs-io)
set_property(TARGET gbs-io PROPERTY EXPORT_NAME gbs::io)

# fmt
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    find_package(fmt REQUIRED)
endif()
# gbs-core
if(NOT TARGET gbs)
    find_package(gbs REQUIRED)
endif ()

# TODO: package magic enum
# TODO: refactor libiges to remove hard-coded incl dirs
target_include_directories(
    gbs-io
    INTERFACE
        $<BUILD_INTERFACE:${GBS_IO_INCLUDE_DIR}>
        $<INSTALL_INTERFACE:include>
        ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}/iges
)

target_link_libraries(
    gbs-io
    INTERFACE
        gbs::core
        iges
)

# link against fmt lib since GCC12 doesn't provide <format> header defined in C++20
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_link_libraries(
        gbs-io
        INTERFACE
            fmt::fmt
    )
endif()

install(FILES ${GBS_IO_HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/gbs/io)
install(FILES ${GBS_IO_3RD_PARTIES} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/gbs/third_parties)

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
