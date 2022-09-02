cmake_minimum_required(VERSION 3.20.0)


# add_definitions(-D_USE_MATH_DEFINES)
find_package(GTest REQUIRED)

include_directories(${CMAKE_SOURCE_DIR})

file(GLOB SRC_LIST "tests_*")
foreach(file_name ${SRC_LIST})
        get_filename_component(exe_name ${file_name} NAME_WE  )
        add_executable(${exe_name} ${file_name})


        target_link_libraries(${exe_name} 
                GTest::gtest
                GTest::gtest_main
                ${VTK_LIBRARIES}
                gbs-render 
        )

        target_compile_options(
                ${exe_name}  PRIVATE 
                $<$<COMPILE_LANGUAGE:CUDA>:
                --extended-lambda 
                # --expt-relaxed-constexpr
                -allow-unsupported-compiler
                >
            )

    if ( MSVC )
        SET_OUT_DIR(${exe_name} ${CMAKE_BINARY_DIR})
    endif ( MSVC )
        
endforeach()