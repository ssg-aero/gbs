# A VTK built with MSVC (conda-forge win-64) exports the MSVC option -utf-8 on VTK::fmt with no
# compiler condition: VTK's ThirdParty/fmt tests the compiler when VTK is built, not when it is
# used. Another compiler then gets it on every file linking VTK (clang's GNU driver reads it as
# `-u tf-8` and warns). Keep it for MSVC only, the condition VTK meant.
#
# Included after find_package(VTK) by the build of gbs, and by GBSConfig.cmake for the projects
# using gbs-render.
if(TARGET VTK::fmt)
    get_target_property(_gbs_vtk_fmt VTK::fmt ALIASED_TARGET)
    if(NOT _gbs_vtk_fmt)
        set(_gbs_vtk_fmt VTK::fmt)
    endif()
    get_target_property(_gbs_vtk_fmt_options ${_gbs_vtk_fmt} INTERFACE_COMPILE_OPTIONS)
    if(_gbs_vtk_fmt_options AND "-utf-8" IN_LIST _gbs_vtk_fmt_options)
        list(REMOVE_ITEM _gbs_vtk_fmt_options "-utf-8")
        list(APPEND _gbs_vtk_fmt_options
            "$<$<OR:$<COMPILE_LANG_AND_ID:C,MSVC>,$<COMPILE_LANG_AND_ID:CXX,MSVC>>:-utf-8>")
        set_target_properties(${_gbs_vtk_fmt} PROPERTIES INTERFACE_COMPILE_OPTIONS "${_gbs_vtk_fmt_options}")
    endif()
    unset(_gbs_vtk_fmt)
    unset(_gbs_vtk_fmt_options)
endif()
