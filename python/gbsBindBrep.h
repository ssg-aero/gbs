#pragma once
// Helpers shared by the bindings of the native BREP kernel (gbsBindBrep.cpp, gbsBindBrepOps.cpp).
#include <pybind11/pybind11.h>

#include <gbs-brep/builders.h>

#include <utility>

namespace gbs_bind
{
    // Builders return std::expected in C++; in Python a failure raises gbs.brep.BRepError, with its
    // code and the entities involved (attributes code and shapes).
    template <typename V>
    V value_or_raise(gbs::brep::BuildResult<V> r)
    {
        return gbs::brep::unwrap(std::move(r));
    }

    // The long operations release the GIL: other Python threads run meanwhile, but none may modify
    // the model they work on.
    using release_gil = pybind11::call_guard<pybind11::gil_scoped_release>;
}
