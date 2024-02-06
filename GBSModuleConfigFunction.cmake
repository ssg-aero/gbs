# Configure GBS modules
# Global list to store all module targets
set(GBS_MODULES "")
# Function to add a C++20 module with dependencies
function(add_cpp20_module module_name)
    # Parse function arguments
    set(oneValueArgs "")
    set(multiValueArgs FILES DEPS DIR)
    cmake_parse_arguments(MODULE "" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

    message(STATUS "Adding module ${module_name}")
    # Add the module
    add_library(${module_name})
    add_library(GBS::${module_name} ALIAS ${module_name})

    # Set module files
    if(MODULE_DIR)
        set(MODULE_FILES_WITH_PATH "")

        foreach(FILE_NAME IN LISTS MODULE_FILES)
            list(APPEND MODULE_FILES_WITH_PATH ${MODULE_DIR}/${FILE_NAME})
        endforeach(FILE_NAME)
    else()
        set(MODULE_FILES_WITH_PATH MODULE_FILES)
    endif()


    # Check for  Clang or GNU version greater than 18/14 for module scan support
    if(
        (CMAKE_CXX_COMPILER_ID STREQUAL "Clang" AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER "18") OR
        (CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER "14")
        )
        target_sources(
            ${module_name} 
            PUBLIC
                FILE_SET 
                    cxx_modules 
                TYPE 
                    CXX_MODULES 
                FILES 
                    ${MODULE_FILES_WITH_PATH}
        )

        set_target_properties(${module_name} PROPERTIES CXX_SCAN_FOR_MODULES on)

    # Check for MSVC
    else()
        target_sources(
            ${module_name}
            PRIVATE
                ${MODULE_FILES_WITH_PATH}
        )
    endif()

    # add_library(${module_name} FILE_SET CXX_MODULES ${MODULE_FILES_WITH_PATH})
    # add_library(GBS::${module_name} ALIAS ${module_name})
    message(STATUS "Composed by files: ${MODULE_FILES_WITH_PATH}")

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
    install(FILES ${MODULE_FILES_WITH_PATH}
        # DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}/gbs/gbs/${module_name}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/gbs/gbs
    )

    # Append to global list
    list(APPEND GBS_MODULES GBS::${module_name})
    set(GBS_MODULES "${GBS_MODULES}" PARENT_SCOPE)
endfunction()