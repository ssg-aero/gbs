cmake_minimum_required(VERSION 3.17.0)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

project(gbs LANGUAGES CXX)

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
        set_property(TARGET ${TARGET_NAME} PROPERTY RUNTIME_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${OUT_DIR} )
        set_property(TARGET ${TARGET_NAME} PROPERTY LIBRARY_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${OUT_DIR} )
        set_property(TARGET ${TARGET_NAME} PROPERTY ARCHIVE_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${OUT_DIR} )
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

message(${LIBRARY_LIB} ${LIBRARY_BIN} ${LIBRARY_INC})
#
# set install dir
#
function(SET_INSTALL_DIR REQUIRED_ARG)
    list(GET ARGV 0 TARGET_NAME)
    if(WIN32)
        if(TARGET ${TARGET_NAME} )
            install(TARGETS ${TARGET_NAME}
                EXPORT  ${TARGET_NAME}Targets
                COMPONENT python
                LIBRARY  DESTINATION ${LIBRARY_LIB}
                RUNTIME  DESTINATION ${LIBRARY_BIN}
                INCLUDES DESTINATION ${LIBRARY_INC}
                ARCHIVE  DESTINATION ${LIBRARY_LIB}
            )
            install(
                EXPORT  ${TARGET_NAME}Targets
                FILE ${TARGET_NAME}.cmake
                NAMESPACE GBS::
                DESTINATION ${LIBRARY_LIB}/cmake/${TARGET_NAME}
            )
        endif(TARGET ${TARGET_NAME} )
    else() # TODO find export dir for linux
            install(TARGETS ${TARGET_NAME}  
            LIBRARY DESTINATION $ENV{LD_RUN_PATH}
            )
    endif()
    # copy headers keeping file structure
    string(TOLOWER ${PROJECT_NAME} inc_sub_dir)
    install(
            DIRECTORY "${CMAKE_SOURCE_DIR}/${inc_sub_dir}"
            DESTINATION "${LIBRARY_INC}"
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

include_directories(${LIBRARY_INC})
link_directories(${LIBRARY_BIN})
link_directories(${LIBRARY_LIB})

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
#
# Boost
# 
# set(BOOST_DIR ${LIBRARY_LIB}/cmake/Boost-1.74.0)
# set(Boost_INCLUDE_DIR ${LIBRARY_INC})
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

SET_INSTALL_DIR(${PROJECT_NAME})

if(${GBS_BUILD_TESTS} OR ${USE_OCCT_UTILS})
    add_subdirectory(gbs-occt)
endif(${GBS_BUILD_TESTS} OR ${USE_OCCT_UTILS})

if(${USE_RENDER})
    add_subdirectory(gbs-render)
endif(${USE_RENDER})

if(${GBS_BUILD_TESTS})
    add_subdirectory(tests)
endif(${GBS_BUILD_TESTS})
if(${USE_PYTHON_BIDING})
    add_subdirectory(python)
endif(${USE_PYTHON_BIDING})