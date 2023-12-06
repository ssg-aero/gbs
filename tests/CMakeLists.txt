cmake_minimum_required(VERSION 3.20.0)

project(gbslib_tests LANGUAGES CXX)

add_definitions(-D_USE_MATH_DEFINES)
# add_definitions(-DTEST_PLOT_ON) # Remvove comments to activate plots

# gslib
include_directories(${CMAKE_SOURCE_DIR})

file(GLOB SRC_LIST "tests_*.cpp")
foreach(file_name ${SRC_LIST})
        get_filename_component(exe_name ${file_name} NAME_WE  )
        add_executable(${exe_name} ${file_name})
#        if(${USE_PROFILE})
#            SET_TARGET_PROPERTIES(${exe_name} PROPERTIES LINK_FLAGS "/PROFILE")
#        endif(${USE_PROFILE})
        
        target_link_libraries(${exe_name} 
        GTest::gtest
        GTest::gtest_main
        NLopt::nlopt
        ${VTK_LIBRARIES}
        gbs-render
        fmt
        iges
        TBB::tbb 
        vecop
        math
        )
        
        install(TARGETS ${exe_name})
        if ( MSVC )
                SET_OUT_DIR(${exe_name} ${CMAKE_BINARY_DIR})
        endif ( MSVC )
        if ( MSVC )
                set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /bigobj")
        endif ( MSVC )
        vtk_module_autoinit(
                TARGETS ${exe_name}
                MODULES ${VTK_LIBRARIES}
                )
        
endforeach()