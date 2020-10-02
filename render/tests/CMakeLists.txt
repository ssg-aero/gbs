project(gbs_render_tests LANGUAGES CXX)
#
# GTest
#
enable_testing()
include_directories(${CONDA_PREFIX}/Library/include)
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