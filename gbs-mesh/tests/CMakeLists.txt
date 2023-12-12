find_package(GTest REQUIRED)

# Set include directories for targets
target_include_directories(gbs-render PUBLIC ${CMAKE_SOURCE_DIR})

foreach(file_name
        tests_halfedgemesh_host.cpp
        tests_halthedgemesh.cu
)
        get_filename_component(exe_name ${file_name} NAME_WE)

        add_executable(${exe_name})
        target_sources(${exe_name} PRIVATE ${file_name})
        target_link_libraries(${exe_name}
                GTest::gtest
                GTest::gtest_main
                ${VTK_LIBRARIES}
                gbs-render
        )

        target_compile_options(
                ${exe_name} PRIVATE
                $<$<COMPILE_LANGUAGE:CUDA>:
                --extended-lambda
                -allow-unsupported-compiler
                >
        )

        set_target_properties(
                ${exe_name} PROPERTIES
                CUDA_SEPARABLE_COMPILATION ON
        )

        if(MSVC)
                SET_OUT_DIR(${exe_name} ${CMAKE_BINARY_DIR})
        endif(MSVC)

        add_test(
                NAME ${exe_name}
                COMMAND ${exe_name}
        )
endforeach()