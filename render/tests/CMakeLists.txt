project(gbs_render_tests LANGUAGES CXX)
#
# GTest
#
enable_testing()
include_directories(${CONDA_PREFIX}/Library/include)
# CODAN ENV
file(TO_CMAKE_PATH $ENV{CONDA_PREFIX} CONDA_PREFIX)
#
# occt
# 
set(OpenCASCADE_DIR ${CONDA_PREFIX}/Library/lib/cmake/opencascade)
find_package(OpenCASCADE REQUIRED)
link_directories(${OpenCASCADE_LIBRARY_DIR})
message(${OpenCASCADE_LIBRARY_DIR})
include_directories(${OpenCASCADE_INCLUDE_DIR})
add_compile_definitions(_SILENCE_CXX17_ITERATOR_BASE_CLASS_DEPRECATION_WARNING) # occt iterators seems to be old
#
# current tested lib
# 
link_directories(${CMAKE_BINARY_DIR})
#
# add sources to test exe
# 
file(GLOB SRC_LIST "tests_*")
add_executable(${PROJECT_NAME} ${SRC_LIST})
target_link_libraries(${PROJECT_NAME} 
                            # GTest::GTest 
                            # GTest::Main 
                            optimized gtest_dll
                            debug gtest_dlld
                            # optimized gtest_main
                            # debug gtest_maind
                            ${OpenCASCADE_LIBRARIES} 
                            gbs_render 
                            nlopt
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