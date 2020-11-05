cmake_minimum_required(VERSION 3.17.0)

project(gbslib_tests LANGUAGES CXX)

# gslib
include_directories(${CMAKE_SOURCE_DIR})

# gbs-occt
include_directories(${CMAKE_SOURCE_DIR}/gbs-occt)

file(GLOB SRC_LIST "tests_*")
add_executable(${PROJECT_NAME} ${SRC_LIST})

target_link_libraries(${PROJECT_NAME} 
                            optimized gtest_dll
                            debug gtest_dlld
                            ${OpenCASCADE_LIBRARIES} 
                            gbs-occt
                            sundials_kinsol
                            nlopt
                            ${VTK_LIBRARIES}
                            gbs-render
                        )

install(TARGETS ${PROJECT_NAME})

vtk_module_autoinit(
        TARGETS ${PROJECT_NAME}
        MODULES ${VTK_LIBRARIES}
)
