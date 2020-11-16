project(gbs-occt_tests LANGUAGES CXX)

file(GLOB SRC_LIST "tests_*")
add_executable(${PROJECT_NAME} ${SRC_LIST})
target_link_libraries(${PROJECT_NAME} 
                            optimized gtest_dll
                            debug gtest_dlld
                            ${OpenCASCADE_LIBRARIES} 
                            gbs-occt 
                    )
install(TARGETS ${PROJECT_NAME})