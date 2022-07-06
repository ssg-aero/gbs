cmake_minimum_required(VERSION 3.20.0)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/)

# Build options
option(USE_CUDA "Enable the CUDA implementation" OFF)
option(BUILD_GBS "Enable the build of gbs-core library" ON)
option(BUILD_GBS_CUDA "Enable the build of gbs-cuda extension" OFF)
option(BUILD_GBS_IO "Enable the build of gbs-io extension" OFF)
option(BUILD_GBS_MESH "Enable the build of gbs-mesh extension" OFF)
option(BUILD_GBS_RENDER "Enable the build of gbs-render extension" OFF)
option(BUILD_GBS_OCCT "Enable the build of gbs-occt extension" OFF)
option(BUILD_PYTHON_BINDINGS "Enable the build of Python bindings" OFF)
option(BUILD_DOC "Enable the build of the C++ API documentation" OFF)
option(BUILD_GBS_TESTS "Enable the build of gbs C++ tests" OFF)

project(gbs LANGUAGES CXX)

# C++ conf
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)

# Set installation directories (CMAKE_INSTALL_INCLUDEDIR, CMAKE_INSTALL_BINDIR, CMAKE_INSTALL_LIBDIR)
include(GNUInstallDirs)

message("-- CMake install prefix: " ${CMAKE_INSTALL_PREFIX})
message("  -> binaries: " ${CMAKE_INSTALL_BINDIR})
message("  -> libs: " ${CMAKE_INSTALL_LIBDIR})
message("  -> includes: " ${CMAKE_INSTALL_INCLUDEDIR})

set(CMAKECONFIG_INSTALL_DIR "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}" CACHE
    STRING "install path for gbs CMake config files")

# gbs-core library
if(${BUILD_GBS})
    add_subdirectory(gbs-core)
endif()

# gbs-cuda extension
if(${BUILD_GBS_CUDA})
    add_subdirectory(gbs-cuda)
endif()

# gbs-io extension
if(${BUILD_GBS_IO})
    add_subdirectory(gbs-io)
endif()

# gbs-mesh extension
if(${BUILD_GBS_MESH})
    add_subdirectory(gbs-mesh)
endif()

# gbs-occt extension
if(${BUILD_GBS_OCCT})
    add_subdirectory(gbs-occt)
endif()

# gbs-render extension
if(${BUILD_GBS_RENDER})
    add_subdirectory(gbs-render)
endif()

# Python bindings
if(${BUILD_PYTHON_BINDINGS})
    add_subdirectory(python)
endif()

# Doxygen documentation
if(${BUILD_DOC})
    add_subdirectory(docs)
endif()

# C++ tests of core and ext libs
if(${BUILD_GBS_TESTS})
    add_subdirectory(tests)
endif()
