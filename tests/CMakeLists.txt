project(gbslib_tests LANGUAGES CXX)

include(GoogleTest)

# add_definitions(-DTEST_PLOT_ON) # Remvove comments to activate plots

# Populate tests list
file(GLOB SRC_LIST "tests_*.cpp")

find_package(libIGES REQUIRED)

foreach(file_name ${SRC_LIST})
        get_filename_component(test_exe_name ${file_name} NAME_WE  )
        add_executable(${test_exe_name} ${file_name})
#        if(${USE_PROFILE})
#            SET_TARGET_PROPERTIES(${exe_name} PROPERTIES LINK_FLAGS "/PROFILE")
#        endif(${USE_PROFILE})
        
        target_link_libraries(${test_exe_name}
                PRIVATE
                        GTest::gtest
                        GTest::gtest_main
                        NLopt::nlopt
                        Eigen3::Eigen
                        ${VTK_LIBRARIES}
                        gbs-render
                        libIGES::iges
                        TBB::tbb
                        ${GBS_MODULES}
        )

        target_include_directories(${test_exe_name}
                PRIVATE
                        ${GTEST_INCLUDE_DIRS}
        )
        
        gtest_discover_tests(${test_exe_name})
        install(TARGETS ${test_exe_name})
        if ( MSVC )
                SET_OUT_DIR(${test_exe_name} ${CMAKE_BINARY_DIR})
        endif ( MSVC )
        if ( MSVC )
                set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /bigobj")
        endif ( MSVC )
        vtk_module_autoinit(
                TARGETS ${test_exe_name}
                MODULES ${VTK_LIBRARIES}
        )
        
endforeach()