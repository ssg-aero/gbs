project(gbs-render_examples LANGUAGES CXX)

# Exemples de rendu/interaction VTK : executables independants, sans framework
# de test (scenes interactives a lancer a la main). Plus de dependance GoogleTest.
file(GLOB EXAMPLE_SRC_LIST "*.cpp")

foreach(example_src ${EXAMPLE_SRC_LIST})
    get_filename_component(example_name ${example_src} NAME_WE)
    add_executable(${example_name} ${example_src})
    target_link_libraries(${example_name}
                                ${OpenCASCADE_LIBRARIES}
                                gbs-render
                                NLopt::nlopt
                                ${VTK_LIBRARIES}
                                TBB::tbb
                        )
    install(TARGETS ${example_name})
    if ( MSVC )
            SET_OUT_DIR(${example_name} ${CMAKE_BINARY_DIR})
    endif ( MSVC )
    vtk_module_autoinit(
            TARGETS ${example_name}
            MODULES ${VTK_LIBRARIES}
    )
endforeach()
