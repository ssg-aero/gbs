project(gbs-occt_tests LANGUAGES CXX)

# tests_*.cpp uniquement (exclut tests_main.cpp.hide, desactive)
file(GLOB SRC_LIST "tests_*.cpp")
add_executable(${PROJECT_NAME} ${SRC_LIST})
target_link_libraries(${PROJECT_NAME}
                        gbs_testing
                        ${OpenCASCADE_LIBRARIES}
                        gbs-occt
                        TBB::tbb
                    )
install(TARGETS ${PROJECT_NAME})
if ( MSVC )
        SET_OUT_DIR(${PROJECT_NAME} ${CMAKE_BINARY_DIR})
endif ( MSVC )
add_test(NAME ${PROJECT_NAME} COMMAND ${PROJECT_NAME})
