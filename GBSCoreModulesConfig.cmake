include("${CMAKE_CURRENT_LIST_DIR}/GBSModuleConfigFunction.cmake")

include(GNUInstallDirs)

find_path(GBS_MODULE_UNITS_PATH NAMES vecop.cppm PATHS gbs/modules ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}/gbs/gbs/modules)

if(GBS_MODULE_UNITS_PATH)
    message(STATUS "Found vecop.cppm in: ${GBS_MODULE_UNITS_PATH}")
else()
    message(STATUS "vecop.cppm not found")
endif()

add_cpp20_module(vecop
    FILES
        "vecop.cppm"
    DEPS
        gbs
    DIR
        ${GBS_MODULE_UNITS_PATH}
)

add_cpp20_module(math
    FILES
        "math.cppm"
    DEPS
        gbs
        GBS::vecop
    DIR
        ${GBS_MODULE_UNITS_PATH}
)

add_cpp20_module(basis_functions
    FILES
        "basisfunctions.cppm"
    DEPS
        gbs
        GBS::vecop
        GBS::math
    DIR
        ${GBS_MODULE_UNITS_PATH}
)

add_cpp20_module(knots_functions
    FILES
        "knotsfunctions.cppm"
    DEPS
        gbs
        GBS::vecop
        GBS::math
        GBS::basis_functions
    DIR
        ${GBS_MODULE_UNITS_PATH}
)