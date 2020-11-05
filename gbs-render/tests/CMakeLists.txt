project(gbs-render_tests LANGUAGES CXX)

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
                            gbs-render 
                            nlopt
                            ${VTK_LIBRARIES}
                    )
install(TARGETS ${PROJECT_NAME})

vtk_module_autoinit(
        TARGETS ${PROJECT_NAME}
        MODULES ${VTK_LIBRARIES}
)