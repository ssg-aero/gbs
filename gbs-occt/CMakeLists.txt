project(gbs-occt LANGUAGES CXX)
set(CMAKE_INCLUDE_CURRENT_DIR ON)
# gbslib
include_directories(${CMAKE_SOURCE_DIR})

file(GLOB OCCT_UTILS_SRC_LIST "*.cpp")
add_library(${PROJECT_NAME} SHARED ${OCCT_UTILS_SRC_LIST})
target_link_libraries(${PROJECT_NAME} 
    gbs
    ${OpenCASCADE_LIBRARIES}
    TBB::tbb    
)

set_target_properties(${PROJECT_NAME} PROPERTIES EXPORT_NAME occt)
install(TARGETS ${PROJECT_NAME} EXPORT GbsTargets
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
INSTALL_HEADERS(${PROJECT_NAME})
if ( MSVC )
        SET_OUT_DIR(${PROJECT_NAME} ${CMAKE_BINARY_DIR})
endif ( MSVC )
# tests
# add current dir for include in tests
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
if(${GBS_BUILD_TESTS})
    add_subdirectory(tests)
endif(${GBS_BUILD_TESTS})

