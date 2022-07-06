cmake_minimum_required(VERSION 3.20.0)

project(test-gbs LANGUAGES CXX)

find_package(GTest REQUIRED)

file(GLOB SRC_LIST "tests_*")
foreach(file_name ${SRC_LIST})
    get_filename_component(exe_name ${file_name} NAME_WE)
    add_executable(${exe_name} ${file_name})

    target_link_libraries(${exe_name} 
            GTest::gtest
            GTest::gtest_main
            gbs-core
            gbs-io
            gbs-render
            gbs-mesh
    )

    install(TARGETS ${exe_name})       
endforeach()