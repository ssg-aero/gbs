project(gbs-occt_tests LANGUAGES CXX)

file(GLOB SRC_LIST "tests_*.cpp")
add_executable(${PROJECT_NAME} ${SRC_LIST})
# libIGES and NLopt: the native BREP IGES export round trip (tests_brep_iges_occt.cpp)
find_package(libIGES REQUIRED)
target_link_libraries(${PROJECT_NAME}
                        gbs_testing
                        ${OpenCASCADE_LIBRARIES}
                        gbs-occt
                        libIGES::iges
                        NLopt::nlopt
                        TBB::tbb
                    )
if ( MSVC )
        SET_OUT_DIR(${PROJECT_NAME} ${CMAKE_BINARY_DIR})
endif ( MSVC )
add_test(NAME ${PROJECT_NAME} COMMAND ${PROJECT_NAME})
