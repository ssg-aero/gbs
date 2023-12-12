find_package(GTest REQUIRED)

include_directories(${CMAKE_SOURCE_DIR})

file(GLOB SOURCE_FILES *.cu)
foreach(file_name ${SOURCE_FILES})
    get_filename_component(exe_name ${file_name} NAME_WE)
    message("exe: " ${exe_name} " file: " ${file_name})
    add_executable(${exe_name} ${file_name})
    target_link_libraries(${exe_name} 
        GTest::gtest
        GTest::gtest_main
    )

    set_target_properties(
        ${exe_name} PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON
    )

    target_compile_options(
        ${exe_name}  PRIVATE 
        $<$<COMPILE_LANGUAGE:CUDA>:
        --extended-lambda 
        -allow-unsupported-compiler
        >
    )

    if (MSVC)
        SET_OUT_DIR(${exe_name} ${CMAKE_BINARY_DIR})
    endif (MSVC)
        
endforeach()