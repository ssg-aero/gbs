#pragma once
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;


#include <gbs/surfaces>
#include <gbs-io/tojson.h>
using namespace gbs;

template <size_t dim>
inline auto add_ext(const char *func_name)
{
    return std::string( func_name ) + std::to_string( dim) + std::string("d");
}

inline const char *build_rep(const auto &cls)
{
    rapidjson::StringBuffer buffer;
    rapidjson::Writer<rapidjson::StringBuffer> writer(buffer);
    rapidjson::Document d;
    d.SetObject();
    auto crv_val = make_json(cls, d.GetAllocator());
    crv_val.Accept(writer);
    return buffer.GetString();
}


template <typename T, size_t dim, bool rational>
inline auto declare_bssurface(py::module_ &m)
{
       
        using Class = typename std::conditional<rational, BSSurfaceRational<T, dim>, BSSurface<T, dim>>::type;
        using ClassBase =  Surface<T,dim>;

        std::string typestr;
        std::string rationalstr;
        if(rational) rationalstr="Rational";
        T x;

        std::string pyclass_name = std::string("BSSurface")+ rationalstr + std::to_string(dim) + "d";// + typestr;
        std::string pyclass_base_name = std::string("BSSurfaceBase")+ rationalstr + std::to_string(dim) + "d";// + typestr;

        py::class_<BSSurfaceGeneral<T, dim, rational>, std::shared_ptr<BSSurfaceGeneral<T, dim, rational>>, ClassBase>(m, pyclass_base_name.c_str());

        return  py::class_<Class,std::shared_ptr<Class>, BSSurfaceGeneral<T, dim, rational>>(m, pyclass_name.c_str())
        .def(py::init<
                        const points_vector<T, dim + rational> &,
                        const std::vector<T> &,
                        const std::vector<T> &,
                        size_t,
                        size_t>(),
                py::arg("poles"), 
                py::arg("knots_flatsU"), 
                py::arg("knots_flatsV"), 
                py::arg("degreeU"),
                py::arg("degreeV")
        )
        .def(py::init<
                        const points_vector<T, dim + rational> &,
                        const std::vector<T> &,
                        const std::vector<T> &,
                        const std::vector<size_t> &,
                        const std::vector<size_t> &,
                        size_t,
                        size_t>(),
                py::arg("poles"), 
                py::arg("knotsU"), 
                py::arg("knotsV"), 
                py::arg("multsU"), 
                py::arg("multsV"), 
                py::arg("degreeU"),
                py::arg("degreeV")
        )
        .def("degreeU", &Class::degreeU,"Surface's U degree")
        .def("degreeV", &Class::degreeV,"Surface's V degree")
        .def("knotsFlatsU", &Class::knotsFlatsU,"Surface's U knots")
        .def("knotsFlatsV", &Class::knotsFlatsV,"Surface's V knots")
        .def("knotsU", &Class::knotsU,"Surface's U knots")
        .def("knotsV", &Class::knotsV,"Surface's V knots")
        .def("multsU", &Class::multsU,"Surface's U knots multiplicities")
        .def("multsV", &Class::multsV,"Surface's V knots multiplicities")
        .def("poles", &Class::poles,"Surface's poles")
        .def("bounds", &Class::bounds,"Returns surface's start stop values")
        .def("insertKnotU", &Class::insertKnotU,"Insert U knot at multiplicity m",py::arg("u"), py::arg("m")=1)
        .def("insertKnotV", &Class::insertKnotV,"Insert V knot at multiplicity m",py::arg("v"), py::arg("m")=1)
        .def("isoU", &Class::isoU,"Isoparametric curve u constant",py::arg("u"))
        .def("isoV", &Class::isoV,"Isoparametric curve v constant",py::arg("v"))
        .def("invertUV",&Class::invertUV,"Invert U, V parametrization")
        .def("trimU",&Class::trimU,"Permanently trim surface between u1 and u2 along V direction",py::arg("u1"),py::arg("u2"))
        .def("trimV",&Class::trimV,"Permanently trim surface between v1 and v2 along U direction",py::arg("v1"),py::arg("v2"))
        .def("__copy__", [](const Class &self)
                { return Class(self); })
        .def(
            py::pickle(
                    [](const Class &srf) {
                            return py::make_tuple(// __getstate__
                                    srf.poles(),
                                    srf.knotsU(),
                                    srf.knotsV(),
                                    srf.multsU(),
                                    srf.multsV(),
                                    srf.degreeU(),
                                    srf.degreeV()
                            );
                    },
                    [](py::tuple t) { // __setstate__
                            if (t.size() != 7)
                                    throw std::runtime_error("Invalid state!");
                            return Class{
                                    t[0].cast<points_vector<T,dim+rational>>(),
                                    t[1].cast<std::vector<T>>(),
                                    t[2].cast<std::vector<T>>(),
                                    t[3].cast<std::vector<size_t>>(),
                                    t[4].cast<std::vector<size_t>>(),
                                    t[5].cast<size_t>(),
                                    t[6].cast<size_t>(),
                            };
                    }
            )
        )
        .def("__repr__", [](const Class &self) { return build_rep( self ); } )
    ;
}

template <typename T, size_t dim>
inline void gbs_bind_surfaces(py::module &m)
{

    py::class_<Surface<T, dim>, std::shared_ptr<Surface<T, dim>> >(m, add_ext<dim>("Surface").c_str())
    .def("value", &Surface<T, dim>::value,"Curve evaluation at given parameter",py::arg("u"),py::arg("v"),py::arg("du") = 0,py::arg("dv") = 0)
    .def("bounds", &Surface<T, dim>::bounds, "Returns surface's start stop values")
    .def("__call__",&Surface<T, dim>::operator(),"Curve evaluation at given parameter",py::arg("u"),py::arg("v"),py::arg("du") = 0,py::arg("dv") = 0)
    ;

    declare_bssurface<T,dim,false>(m);
    declare_bssurface<T,dim,true>(m)
    .def(
            py::init<
                const points_vector<T, dim> &, 
                const std::vector<T> &, 
                const std::vector<T> &, 
                const std::vector<T> &, 
                size_t, 
                size_t
            >(), 
            py::arg("poles"), 
            py::arg("weights"), 
            py::arg("knots_flatsU"), 
            py::arg("knots_flatsV"), 
            py::arg("degreeU"),
            py::arg("degreeV")
    )
    .def(
        py::init<const points_vector<T, dim> &, const std::vector<T> &, const std::vector<T> &, const std::vector<T> &, const std::vector<size_t> &, const std::vector<size_t> &, size_t, size_t>(), 
        py::arg("poles"), 
        py::arg("weights"), 
        py::arg("knotsU"), 
        py::arg("knotsV"), 
        py::arg("multsU"), 
        py::arg("multsV"), 
        py::arg("degreeU"),
        py::arg("degreeV")
    )
    .def("polesProjected",&BSSurfaceRational<T, dim>::polesProjected)
    .def("weights",&BSSurfaceRational<T, dim>::weights)
    ;

}

