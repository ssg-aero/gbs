project(gbs-render_tests LANGUAGES CXX)

file(GLOB SRC_LIST "tests_*")
add_executable(${PROJECT_NAME} ${SRC_LIST})
target_link_libraries(${PROJECT_NAME} 
                            GTest::gtest
                            GTest::gtest_main
                            ${OpenCASCADE_LIBRARIES} 
                            gbs-render 
                            NLopt::nlopt
                            ${VTK_LIBRARIES}
                    )
install(TARGETS ${PROJECT_NAME})
if ( MSVC )
        SET_OUT_DIR(${PROJECT_NAME} ${CMAKE_BINARY_DIR})
endif ( MSVC )
vtk_module_autoinit(
        TARGETS ${PROJECT_NAME}
        MODULES ${VTK_LIBRARIES}
)