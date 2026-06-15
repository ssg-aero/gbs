# Test host (CPU pur) : construit des que GBS_BUILD_TESTS est actif, sans CUDA,
# pour qu'il tourne en CI. Le test CUDA reste opt-in (GBS_USE_CUDA) car il exige
# un GPU absent des runners GitHub.
set(GBS_MESH_TEST_FILES tests_halfedgemesh_host.cpp)
if(GBS_USE_CUDA)
        list(APPEND GBS_MESH_TEST_FILES tests_halthedgemesh.cu)
endif()

foreach(file_name ${GBS_MESH_TEST_FILES})
        get_filename_component(exe_name ${file_name} NAME_WE)

        add_executable(${exe_name})
        target_sources(${exe_name} PRIVATE ${file_name})
        target_include_directories(${exe_name} PRIVATE ${CMAKE_SOURCE_DIR})
        target_link_libraries(${exe_name}
                gbs_testing
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