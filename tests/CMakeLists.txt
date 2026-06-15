project(gbslib_tests LANGUAGES CXX)

# doctest fournit son module de decouverte de tests a cote de sa config.
include(${doctest_DIR}/doctest.cmake)

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
                        gbs_testing
                        NLopt::nlopt
                        Eigen3::Eigen
                        ${VTK_LIBRARIES}
                        gbs-render
                        libIGES::iges
                        TBB::tbb
                        ${GBS_MODULES}
        )

        doctest_discover_tests(${test_exe_name})
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

# Source-based coverage target (only meaningful when configured with
# -DGBS_COVERAGE=ON). Runs the priority-perimeter test executables and
# aggregates an llvm-cov report under ${CMAKE_BINARY_DIR}/coverage/.
# Reproduce from scratch with scripts/coverage.sh.
if(GBS_COVERAGE)
        set(GBS_COVERAGE_TESTS
                tests_basisfunctions
                tests_bscurve
                tests_bssurf
                tests_curve_derivatives
                tests_surface_derivatives
                tests_curves
                tests_surfaces
                tests_bsinterp
                tests_knotsfunctions
                tests_bssbuild
                tests_bsctools
        )
        add_custom_target(coverage
                COMMAND ${CMAKE_COMMAND} -E env bash
                        ${CMAKE_SOURCE_DIR}/scripts/coverage_report.sh
                        ${CMAKE_BINARY_DIR} ${GBS_COVERAGE_TESTS}
                DEPENDS ${GBS_COVERAGE_TESTS}
                WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
                COMMENT "Running instrumented tests and aggregating llvm-cov report"
                USES_TERMINAL
        )
endif()
