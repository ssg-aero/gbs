cmake_minimum_required(VERSION 3.17.0)

project(gbslib_tests LANGUAGES CXX)
#
# c++ conf
#
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)
set(CMAKE_INCLUDE_CURRENT_DIR ON)
#
# gslib
#
message(${CMAKE_CURRENT_SOURCE_DIR}/../)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../)
#
# CODAN ENV
#
file(TO_CMAKE_PATH $ENV{CONDA_PREFIX} CONDA_PREFIX)
include_directories(${CONDA_PREFIX}/Library/include)
link_directories(${CONDA_PREFIX}/Library/lib)
link_directories(${CONDA_PREFIX}/Library/bin)
#
# GTest
#
enable_testing()
include_directories(${CONDA_PREFIX}/Library/include)
#
# eigen
#
include_directories(${CONDA_PREFIX}/Library/include/eigen3)
#
# occt
# 
set(OpenCASCADE_DIR ${CONDA_PREFIX}/Library/lib/cmake/opencascade)
find_package(OpenCASCADE REQUIRED)
link_directories(${OpenCASCADE_LIBRARY_DIR})
include_directories(${OpenCASCADE_INCLUDE_DIR})
add_compile_definitions(_SILENCE_CXX17_ITERATOR_BASE_CLASS_DEPRECATION_WARNING) # occt interators seems to be old
#
# occt-utils
#
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../occt-utils)
# 
# VTK
# 
set(VTK_DIR  ${CONDA_PREFIX}/Library/lib/cmake/vtk-9.0)
include_directories(${CONDA_PREFIX}/Library/include/vtk-9.0)
find_package(VTK COMPONENTS 
CommonCore
CommonColor
CommonDataModel
FiltersSources
InteractionStyle
RenderingCore
RenderingFreeType
RenderingOpenGL2
)

file(GLOB SRC_LIST "tests_*")
message(${SRC_LIST})
add_executable(${PROJECT_NAME} ${SRC_LIST})

target_link_libraries(${PROJECT_NAME} 
                        #     GTest::GTest 
                        #     GTest::Main 
                        optimized gtest_dll
                        debug gtest_dlld
                        # optimized gtest
                        # debug gtestd
                        # optimized gtest_main
                        # debug gtest_maind
                            ${OpenCASCADE_LIBRARIES} 
                            occt-utils
                            sundials_kinsol
                            nlopt
                            ${VTK_LIBRARIES}
                            gbs_render
                    )

#
# set output dirs (vs is a pain)
#
if ( MSVC )
    foreach( OUTPUTCONFIG ${CMAKE_CONFIGURATION_TYPES} )
        string( TOUPPER ${OUTPUTCONFIG} OUTPUTCONFIG )
        set_property(TARGET ${PROJECT_NAME} PROPERTY RUNTIME_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${CMAKE_BINARY_DIR} )
        set_property(TARGET ${PROJECT_NAME} PROPERTY LIBRARY_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${CMAKE_BINARY_DIR} )
        set_property(TARGET ${PROJECT_NAME} PROPERTY ARCHIVE_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${CMAKE_BINARY_DIR} )
    endforeach( OUTPUTCONFIG CMAKE_CONFIGURATION_TYPES )
endif ( MSVC )

vtk_module_autoinit(
        TARGETS ${PROJECT_NAME}
        MODULES ${VTK_LIBRARIES}
)
