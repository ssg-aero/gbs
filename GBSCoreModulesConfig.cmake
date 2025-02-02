include("${CMAKE_CURRENT_LIST_DIR}/GBSModuleConfigFunction.cmake")

include(GNUInstallDirs)

find_path(MODULE_FILES_PATH NAMES vecop.ixx PATHS gbs ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}/gbs/gbs)

if(MODULE_FILES_PATH)
    message(STATUS "Found vecop.ixx in: ${MODULE_FILES_PATH}")
else()
    message(STATUS "vecop.ixx not found")
endif()

add_cpp20_module(vecop
    FILES
        "vecop.ixx"
    DIR
        ${MODULE_FILES_PATH}
)

add_cpp20_module(math
    FILES
        "math.ixx"
    DEPS
        GBS::vecop
    DIR
        ${MODULE_FILES_PATH}
)

add_cpp20_module(basis_functions
    FILES
        "basisfunctions.ixx"
    DEPS
        GBS::vecop
        GBS::math
    DIR
        ${MODULE_FILES_PATH}
)

add_cpp20_module(knots_functions
    FILES
        "knotsfunctions.ixx"
    DEPS
        GBS::vecop
        GBS::math
        GBS::basis_functions
    DIR
        ${MODULE_FILES_PATH}
)