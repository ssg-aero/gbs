cmake_minimum_required(VERSION 3.20.0)

project(gbslib_tests LANGUAGES CXX)

add_definitions(-D_USE_MATH_DEFINES)


# gslib
include_directories(${CMAKE_SOURCE_DIR})

# gbs-occt
include_directories(${CMAKE_SOURCE_DIR}/gbs-occt)

file(GLOB SRC_LIST "tests_*")
foreach(file_name ${SRC_LIST})
        get_filename_component(exe_name ${file_name} NAME_WE  )
        add_executable(${exe_name} ${file_name})

        target_link_libraries(${exe_name} 
        GTest::gtest
        GTest::gtest_main
        ${OpenCASCADE_LIBRARIES} 
        gbs-occt
        NLopt::nlopt
        ${VTK_LIBRARIES}
        gbs-render
        fmt
        iges
        )

        install(TARGETS ${exe_name})
        if ( MSVC )
        SET_OUT_DIR(${exe_name} ${CMAKE_BINARY_DIR})
        endif ( MSVC )
        vtk_module_autoinit(
                TARGETS ${exe_name}
                MODULES ${VTK_LIBRARIES}
                )
        
endforeach()