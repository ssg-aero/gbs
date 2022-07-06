cmake_minimum_required(VERSION 3.17.0)

project(gbs-io LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)

set(GBS_IO_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(GBS_IO_INCLUDE_3RD_PARTIES ${CMAKE_CURRENT_SOURCE_DIR}/third_parties)

set(GBS_IO_HEADERS
    ${GBS_IO_INCLUDE_DIR}/gbs/io/fromjson.h
    ${GBS_IO_INCLUDE_DIR}/gbs/io/fromtext.h
    ${GBS_IO_INCLUDE_DIR}/gbs/io/iges.h
    ${GBS_IO_INCLUDE_DIR}/gbs/io/print.h
    ${GBS_IO_INCLUDE_DIR}/gbs/io/tojson.h
)
add_library(gbs-io INTERFACE ${GBS_IO_HEADERS})
add_library(gbs::io ALIAS gbs-io)

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
        ${GBS_IO_INCLUDE_3RD_PARTIES}
        ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}/iges
        gbs::core
)
target_link_libraries(
    gbs-io
    INTERFACE
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
