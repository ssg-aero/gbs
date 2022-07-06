project(test-gbs-occt LANGUAGES CXX)

find_package(GTest REQUIRED)

file(GLOB SRC_LIST "tests_*")

add_executable(${PROJECT_NAME} ${SRC_LIST})

target_link_libraries(${PROJECT_NAME} 
                        GTest::gtest
                        GTest::gtest_main
                        gbs::occt
                    )

install(TARGETS ${PROJECT_NAME})
