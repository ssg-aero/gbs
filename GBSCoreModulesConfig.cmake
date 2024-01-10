include("${CMAKE_CURRENT_LIST_DIR}/GBSModuleConfigFunction.cmake")

add_cpp20_module(vecop
    FILES
        "gbs/vecop.ixx"
)

add_cpp20_module(math
    FILES
        "gbs/math.ixx"
    DEPS
        GBS::vecop
)

add_cpp20_module(basis_functions
    FILES
        "gbs/basisfunctions.ixx"
    DEPS
        GBS::vecop
        GBS::math
)

add_cpp20_module(knots_functions
    FILES
        "gbs/knotsfunctions.ixx"
    DEPS
        GBS::vecop
        GBS::math
        GBS::basis_functions
)