cmake_minimum_required(VERSION 3.28) # For Modules
cmake_minimum_required(VERSION 3.28) # For Modules

project(gbs LANGUAGES C CXX)

# c++ conf
set(CMAKE_CXX_STANDARD 20)       

set(CMAKE_CXX_STANDARD_REQUIRED ON)
# if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libstdc++")
    # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
    # include_directories(/usr/include/clang/18)
# endif()
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --fmodules-ts")
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated-declarations")
endif()
set(CMAKE_CXX_EXTENSIONS ON)
set(CMAKE_INCLUDE_CURRENT_DIR ON)
# set(CMAKE_EXPORT_COMPILE_COMMANDS ON) # Usefull for debug
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

function(INSTALL_HEADERS DIR_PATH)
    # copy headers keeping file structure
    string(TOLOWER ${DIR_PATH} inc_sub_dir)
    install(
            DIRECTORY "${CMAKE_SOURCE_DIR}/${inc_sub_dir}"
            DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/gbs"
            FILES_MATCHING 
            PATTERN "*.h"
            PATTERN "*.hpp"
            PATTERN "*.ixx"
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
function(SET_OUT_DIR TARGET_NAME OUT_DIR)
    foreach(OUTPUTCONFIG IN LISTS CMAKE_CONFIGURATION_TYPES)
        string(TOUPPER ${OUTPUTCONFIG} OUTPUTCONFIG)
        set_property(TARGET ${TARGET_NAME} PROPERTY RUNTIME_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${OUT_DIR})
        set_property(TARGET ${TARGET_NAME} PROPERTY LIBRARY_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${OUT_DIR})
        set_property(TARGET ${TARGET_NAME} PROPERTY ARCHIVE_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${OUT_DIR})
    endforeach()
endfunction()

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/inc)

# Boost
find_package(Boost REQUIRED)

# NLOPT
find_package(NLopt REQUIRED)

# Eigen
find_package (Eigen3 REQUIRED NO_MODULE)
message(${EIGEN3_INCLUDE_DIR})
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

if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    # ill detected
    # https://stackoverflow.com/questions/71565682/returning-a-rapidjsongenericvalue-from-a-function
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "17.0")
        message(FATAL_ERROR "Insufficient clang version")
    endif()
    add_definitions(-DRAPIDJSON_HAS_CXX11_RVALUE_REFS)
    set(USING_CLANG ON)
endif()


INSTALL_HEADERS(${PROJECT_NAME})
INSTALL_HEADERS("inc")

#Define defaults
if(NOT DEFINED USE_PYTHON_BINDINGS)
    set(USE_PYTHON_BINDINGS ON)
    message(STATUS "USE_PYTHON_BINDINGS default set to " ${USE_PYTHON_BINDINGS})
endif()
if(NOT DEFINED USE_RENDER)
    set(USE_RENDER OFF)
    message(STATUS "USE_RENDER default set to " ${USE_RENDER})
endif()
if(NOT DEFINED USE_CUDA)
    set(USE_CUDA OFF)
    message(STATUS "USE_CUDA default set to " ${USE_CUDA})
endif()
if(NOT DEFINED GBS_BUILD_TESTS)
    set(GBS_BUILD_TESTS ON)
    message(STATUS "GBS_BUILD_TESTS set to " ${GBS_BUILD_TESTS})
endif()
if(NOT DEFINED GBS_USE_MODULES)
    message(STATUS "GBS c++20 modules not used.")
    set(GBS_USE_MODULES OFF)
    endif()
    
    # Include the separate script for module configuration
    # include(GBSModuleConfigFunction.cmake)
if(${GBS_USE_MODULES})

    add_definitions(-DGBS_USE_MODULES)
    set(CXX_SCAN_FOR_MODULES ON)   
    include(GBSCoreModulesConfig.cmake)

    install( FILES GBSModuleConfigFunction.cmake GBSCoreModulesConfig.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/GBSCoreModules
    )


    # link_directories(${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}/gbs)
    # Generate and install the GBSConfig.cmake file
    # install(EXPORT GBSModules
    #     FILE GBSConfig.cmake
    #     NAMESPACE GBS::
    #     DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/GBS
    # )

    message(STATUS "GBS_MODULES:")
    foreach(module IN LISTS GBS_MODULES)
        message(STATUS " * ${module}")
    endforeach()

endif()

# Tests
if(${GBS_BUILD_TESTS})
    message(STATUS "Configurig tests.")

    # Google Tests #
    if(
        DEFINED ENV{CONDA_BUILD} 
        OR DEFINED ENV{PREFIX} 
        OR DEFINED ENV{BUILD_PREFIX} 
        # OR DEFINED USING_CLANG # Clang issue with google test compilation
    )
        message(STATUS "Detected Conda build environment.")
        find_package(GTest REQUIRED)
        message(STATUS "gtest include dirs: " ${GTEST_INCLUDE_DIRS})
        message(STATUS "gtest lib dirs: " ${GTEST_LIBRARIES})

    else() # Local source code is required to build safly in debug mode
        message(STATUS "Building outside of Conda environment.")
        set(GTEST_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/googletest/include)
        include_directories(${GTEST_INCLUDE_DIRS})
        add_subdirectory(googletest)

        if(NOT TARGET GTest::gtest)
            add_library(GTest::gtest ALIAS gtest)
        endif()
        if(NOT TARGET GTest::gtest_main)
            add_library(GTest::gtest_main ALIAS gtest_main)
        endif()
    endif()


    add_subdirectory(tests)
    if(${GBS_BUILD_TESTS_3RDLIBS})
        add_subdirectory(tests_3rdlibs)
    endif(${GBS_BUILD_TESTS_3RDLIBS})

endif()



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
if(USE_RENDER OR GBS_BUILD_TESTS OR USE_PYTHON_BINDINGS)
    set(USE_RENDER ON)
endif()
if(${USE_RENDER}) 
    message("Configure rendering lib " ${USE_RENDER})
    add_subdirectory(gbs-render)
endif(${USE_RENDER})

if(${USE_CUDA})
    message("Configure CUDA implemetation " ${USE_CUDA})
    add_subdirectory(gbs-cuda)
endif(${USE_CUDA})

add_subdirectory(gbs-io)
add_subdirectory(gbs-mesh)
add_subdirectory(tools)

# occt
if(${USE_OCCT_UTILS})
    message("Configure occt bridge " ${USE_OCCT_UTILS})
    find_package(OpenCASCADE REQUIRED)
    link_directories(${OpenCASCADE_LIBRARY_DIR})
    include_directories(${OpenCASCADE_INCLUDE_DIR})
    add_compile_definitions(_SILENCE_CXX17_ITERATOR_BASE_CLASS_DEPRECATION_WARNING) # occt interators seems to be old
    add_subdirectory(gbs-occt)
endif(${USE_OCCT_UTILS})

# Python bindings
if(${USE_PYTHON_BINDINGS})
    message("Configure python binding " ${USE_PYTHON_BINDINGS})
    add_subdirectory(python)
endif(${USE_PYTHON_BINDINGS})

if(${BUILD_DOC})
    add_subdirectory(docs)
endif(${BUILD_DOC})

# add_executable(main main.cpp)
# SET_TARGET_PROPERTIES(main PROPERTIES LINK_FLAGS "/PROFILE")