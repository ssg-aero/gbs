project(occt-utils_tests LANGUAGES CXX)
#
# add sources to test exe
# 
file(GLOB SRC_LIST "tests_*")
add_executable(${PROJECT_NAME} ${SRC_LIST})
target_link_libraries(${PROJECT_NAME} 
                            optimized gtest_dll
                            debug gtest_dlld
                            ${OpenCASCADE_LIBRARIES} 
                            occt-utils 
                    )
if ( MSVC )
    SET_OUT_DIR(${PROJECT_NAME} ${CMAKE_BINARY_DIR})
endif ( MSVC )