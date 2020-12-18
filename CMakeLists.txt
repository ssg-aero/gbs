cmake_minimum_required(VERSION 3.17.0)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/)

project(gbs LANGUAGES CXX)

# c++ conf
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)
set(CMAKE_INCLUDE_CURRENT_DIR ON)

# Set installation directories (CMAKE_INSTALL_INCLUDEDIR, CMAKE_INSTALL_BINDIR, CMAKE_INSTALL_LIBDIR)
include(GNUInstallDirs)

message("-- CMake install prefix: " ${CMAKE_INSTALL_PREFIX})
message("  -> binaries: " ${CMAKE_INSTALL_BINDIR})
message("  -> libs: " ${CMAKE_INSTALL_LIBDIR})
message("  -> includes: " ${CMAKE_INSTALL_INCLUDEDIR})

function(INSTALL_HEADERS REQUIRED_ARG)
    # copy headers keeping file structure
    string(TOLOWER ${PROJECT_NAME} inc_sub_dir)
    install(
            DIRECTORY "${CMAKE_SOURCE_DIR}/${inc_sub_dir}"
            DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
            FILES_MATCHING PATTERN "*.h"
            PATTERN "cmake" EXCLUDE
            PATTERN "tests" EXCLUDE
            PATTERN "build*" EXCLUDE
            PATTERN "python" EXCLUDE
            PATTERN "xml" EXCLUDE
            PATTERN "doc" EXCLUDE
            PATTERN "recipe" EXCLUDE
            PATTERN ".*" EXCLUDE
    )
endfunction()

#
# set output dirs (vs is a pain)
#
function(SET_OUT_DIR REQUIRED_ARG)
    list(GET ARGV 0 TARGET_NAME)
    list(GET ARGV 1 OUT_DIR)
    message(${TARGET_NAME})
    message(${OUT_DIR})
    foreach( OUTPUTCONFIG ${CMAKE_CONFIGURATION_TYPES} )
        string( TOUPPER ${OUTPUTCONFIG} OUTPUTCONFIG )
        set_property(TARGET ${TARGET_NAME} PROPERTY RUNTIME_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${OUT_DIR} )
        set_property(TARGET ${TARGET_NAME} PROPERTY LIBRARY_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${OUT_DIR} )
        set_property(TARGET ${TARGET_NAME} PROPERTY ARCHIVE_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${OUT_DIR} )
    endforeach( OUTPUTCONFIG CMAKE_CONFIGURATION_TYPES )
endfunction()

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# Boost
find_package(Boost REQUIRED)

# Eigen
find_package (Eigen3 3.3 REQUIRED NO_MODULE)
include_directories(${EIGEN3_INCLUDE_DIR})

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
            RenderingAnnotation)

# occt
find_package(OpenCASCADE REQUIRED)
link_directories(${OpenCASCADE_LIBRARY_DIR})
include_directories(${OpenCASCADE_INCLUDE_DIR})
add_compile_definitions(_SILENCE_CXX17_ITERATOR_BASE_CLASS_DEPRECATION_WARNING) # occt interators seems to be old

INSTALL_HEADERS(${PROJECT_NAME})

# Tests
if(${GBS_BUILD_TESTS})
    # Google Tests #
    find_package(GTest REQUIRED)
    include_directories(${GTEST_INCLUDE_DIRS})
    link_directories(${GTEST_LIBRARIES})

    message(STATUS "gtest include dirs: " ${GTEST_INCLUDE_DIRS})
    message(STATUS "gtest lib dirs: " ${GTEST_LIBRARIES})
    
    add_subdirectory(tests)
endif(${GBS_BUILD_TESTS})

# OCCT utils
if(${GBS_BUILD_TESTS} OR ${USE_OCCT_UTILS})
    add_subdirectory(gbs-occt)
endif(${GBS_BUILD_TESTS} OR ${USE_OCCT_UTILS})

# Rendering module
if(${USE_RENDER} OR ${GBS_BUILD_TESTS})
    add_subdirectory(gbs-render)
endif(${USE_RENDER} OR ${GBS_BUILD_TESTS})

# Python bindings
if(${USE_PYTHON_BINDINGS})
    add_subdirectory(python)
endif(${USE_PYTHON_BINDINGS})

if(${BUILD_DOC})
    add_subdirectory(docs)
endif(${BUILD_DOC})