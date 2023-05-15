cmake_minimum_required(VERSION 3.25.0)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/)

project(gbs LANGUAGES C CXX)

# c++ conf
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)
set(CMAKE_INCLUDE_CURRENT_DIR ON)
if ( MSVC )
    # set(CMAKE_CXX_FLAGS     "${CMAKE_CXX_FLAGS} /fp:fast /arch:AVX512")
    set(CMAKE_CXX_FLAGS     "${CMAKE_CXX_FLAGS} /fp:fast /arch:AVX2")
endif( MSVC )

# Set installation directories (CMAKE_INSTALL_INCLUDEDIR, CMAKE_INSTALL_BINDIR, CMAKE_INSTALL_LIBDIR)
include(GNUInstallDirs)

message("-- CMake install prefix: " ${CMAKE_INSTALL_PREFIX})
message("  -> binaries: " ${CMAKE_INSTALL_BINDIR})
message("  -> libs: " ${CMAKE_INSTALL_LIBDIR})
message("  -> includes: " ${CMAKE_INSTALL_INCLUDEDIR})
link_directories(${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})
link_directories(${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR})
include_directories(${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR})

function(INSTALL_HEADERS REQUIRED_ARG)
    # copy headers keeping file structure
    string(TOLOWER ${REQUIRED_ARG} inc_sub_dir)
    install(
            DIRECTORY "${CMAKE_SOURCE_DIR}/${inc_sub_dir}"
            DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/gbs"
            FILES_MATCHING 
            PATTERN "*.h"
            PATTERN "*.hpp"
            PATTERN "curves"
            PATTERN "surfaces"
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
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/inc)

# Boost
find_package(Boost REQUIRED)

# NLOPT
find_package(NLopt REQUIRED)

# Eigen
find_package (Eigen3 3.3 REQUIRED NO_MODULE)
include_directories(${EIGEN3_INCLUDE_DIR})
# TBB
if(UNIX)
    find_package(TBB REQUIRED)
endif()

# libIGES
include_directories(${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}/iges/)

# Python Dev
find_package (Python3 COMPONENTS Interpreter Development)

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
    PythonInterpreter
    WrappingPythonCore
)



INSTALL_HEADERS(${PROJECT_NAME})
INSTALL_HEADERS("inc")

# Tests
if(${GBS_BUILD_TESTS})
    # Google Tests #
    find_package(GTest REQUIRED)
    include_directories(${GTEST_INCLUDE_DIRS})
    link_directories(${GTEST_LIBRARIES})

    message(STATUS "gtest include dirs: " ${GTEST_INCLUDE_DIRS})
    message(STATUS "gtest lib dirs: " ${GTEST_LIBRARIES})
    
    add_subdirectory(tests)
    if(${GBS_BUILD_TESTS_3RDLIBS})
        add_subdirectory(tests_3rdlibs)
    endif(${GBS_BUILD_TESTS_3RDLIBS})
endif(${GBS_BUILD_TESTS})


if(${USE_CUDA})
    find_package(CUDA REQUIRED)
    include_directories(${CUDA_INCLUDE_DIRS})
    enable_language(CUDA)
    set(CMAKE_CUDA_STANDARD 20)
    if(%{CUDA_ARCHITECTURES})
        set(CMAKE_CUDA_ARCHITECTURES %{CUDA_ARCHITECTURES})
    else()
        set(CMAKE_CUDA_ARCHITECTURES 61) # Max 61 for SB2
    endif()
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)
endif()

# Rendering module
if(${USE_RENDER} OR ${GBS_BUILD_TESTS} OR ${USE_PYTHON_BINDINGS}) 
    add_subdirectory(gbs-render)
endif(${USE_RENDER} OR ${GBS_BUILD_TESTS} OR ${USE_PYTHON_BINDINGS})

if(${USE_CUDA})
    add_subdirectory(gbs-cuda)
endif(${USE_CUDA})

add_subdirectory(gbs-io)
add_subdirectory(gbs-mesh)
add_subdirectory(tools)

# occt
if(${USE_OCCT_UTILS})
    find_package(OpenCASCADE REQUIRED)
    link_directories(${OpenCASCADE_LIBRARY_DIR})
    include_directories(${OpenCASCADE_INCLUDE_DIR})
    add_compile_definitions(_SILENCE_CXX17_ITERATOR_BASE_CLASS_DEPRECATION_WARNING) # occt interators seems to be old
    add_subdirectory(gbs-occt)
endif(${USE_OCCT_UTILS})

# # OCCT utils
# if(${GBS_BUILD_TESTS} OR ${USE_OCCT_UTILS})
#     add_subdirectory(gbs-occt)
# endif(${GBS_BUILD_TESTS} OR ${USE_OCCT_UTILS})

# Python bindings
if(${USE_PYTHON_BINDINGS})
    add_subdirectory(python)
endif(${USE_PYTHON_BINDINGS})

if(${BUILD_DOC})
    add_subdirectory(docs)
endif(${BUILD_DOC})

# add_executable(main main.cpp)
# SET_TARGET_PROPERTIES(main PROPERTIES LINK_FLAGS "/PROFILE")