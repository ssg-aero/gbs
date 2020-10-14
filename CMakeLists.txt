cmake_minimum_required(VERSION 3.17.0)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

project(gbslib LANGUAGES CXX)

#
# c++ conf
#
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)
set(CMAKE_INCLUDE_CURRENT_DIR ON)
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
        set_property(TARGET ${PROJECT_NAME} PROPERTY RUNTIME_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${CMAKE_BINARY_DIR} )
        set_property(TARGET ${PROJECT_NAME} PROPERTY LIBRARY_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${CMAKE_BINARY_DIR} )
        set_property(TARGET ${PROJECT_NAME} PROPERTY ARCHIVE_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${CMAKE_BINARY_DIR} )
    endforeach( OUTPUTCONFIG CMAKE_CONFIGURATION_TYPES )
endfunction()
#
# CODAN ENV
#
if(EXISTS $ENV{LIBRARY_PREFIX})
    file(TO_CMAKE_PATH $ENV{LIBRARY_BIN} LIBRARY_BIN)
    file(TO_CMAKE_PATH $ENV{LIBRARY_INC} LIBRARY_INC)
    file(TO_CMAKE_PATH $ENV{LIBRARY_LIB} LIBRARY_LIB)
elseif(EXISTS $ENV{CONDA_PREFIX})
    file(TO_CMAKE_PATH $ENV{CONDA_PREFIX} CONDA_PREFIX)
    file(TO_CMAKE_PATH ${CONDA_PREFIX}/Library/bin LIBRARY_BIN)
    file(TO_CMAKE_PATH ${CONDA_PREFIX}/Library/include LIBRARY_INC)
    file(TO_CMAKE_PATH ${CONDA_PREFIX}/Library/lib LIBRARY_LIB)
else()
    message(FATAL_ERROR "A conda environement is required")
endif()

include_directories(${LIBRARY_INC})
link_directories(${LIBRARY_BIN})
link_directories(${LIBRARY_LIB})

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
#
# Boost
# 
find_package(Boost REQUIRED)
#
# eigen
#
include_directories(${LIBRARY_INC}/eigen3)
# 
# VTK
# 
set(VTK_DIR  ${LIBRARY_LIB}/cmake/vtk-9.0)
include_directories(${LIBRARY_INC}/vtk-9.0)
find_package(VTK COMPONENTS
            CommonCore
            CommonColor
            CommonDataModel
            FiltersSources
            InteractionStyle
            RenderingCore
            RenderingFreeType
            RenderingOpenGL2)
#
# occt
# 
set(OpenCASCADE_DIR ${CONDA_PREFIX}/Library/lib/cmake/opencascade)
find_package(OpenCASCADE )
link_directories(${OpenCASCADE_LIBRARY_DIR})
include_directories(${OpenCASCADE_INCLUDE_DIR})
add_compile_definitions(_SILENCE_CXX17_ITERATOR_BASE_CLASS_DEPRECATION_WARNING) # occt interators seems to be old
#
# current tested lib
# 
link_directories(${CMAKE_BINARY_DIR})
############
if(${GBS_BUILD_TESTS} OR ${USE_OCCT_UTILS})
    add_subdirectory(occt-utils)
endif(${GBS_BUILD_TESTS} OR ${USE_OCCT_UTILS})

if(${USE_RENDER})
    add_subdirectory(render)
endif(${USE_RENDER})

if(${GBS_BUILD_TESTS})
    add_subdirectory(tests)
endif(${GBS_BUILD_TESTS})
if(${USE_PYTHON_BIDING})
    add_subdirectory(python)
endif(${USE_PYTHON_BIDING})