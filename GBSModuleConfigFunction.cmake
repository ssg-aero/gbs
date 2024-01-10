# Configure GBS modules
# Global list to store all module targets
set(GBS_MODULES "")
# Function to add a C++20 module with dependencies
function(add_cpp20_module module_name)
    # Parse function arguments
    set(oneValueArgs "")
    set(multiValueArgs FILES DEPS)
    cmake_parse_arguments(MODULE "" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

    # Add the module
    add_library(${module_name})
    add_library(GBS::${module_name} ALIAS ${module_name})
    message(STATUS "Adding module ${module_name}")
    message(STATUS "Composed by files: ${MODULE_FILES}")

    # Set module files
    target_sources(${module_name}
        # PUBLIC
        #     FILE_SET CXX_MODULES FILES ${MODULE_FILES}
        PRIVATE
            ${MODULE_FILES}
    )

    # Link dependencies if any
    if(MODULE_DEPS)
        target_link_libraries(${module_name} PRIVATE ${MODULE_DEPS})
    endif()

    # add to install
    install(TARGETS ${module_name} 
        # EXPORT GBSModules
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}/gbs
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}/gbs
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/gbs/gbs
        # RENAME "gbs_${module_name}.lib"
    )
    # Install module interface files ${CMAKE_INSTALL_PREFIX}/
    install(FILES ${MODULE_FILES}
        # DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}/gbs/gbs/${module_name}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/gbs/gbs
    )

    # Append to global list
    list(APPEND GBS_MODULES GBS::${module_name})
    set(GBS_MODULES "${GBS_MODULES}" PARENT_SCOPE)
endfunction()