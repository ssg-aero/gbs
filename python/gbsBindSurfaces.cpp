#include "gbsBindSurfaces.h"

void gbs_bind_surfaces(py::module &m)
{
        gbs_bind_surfaces<double,1>(m);
        gbs_bind_surfaces<double,2>(m);
        gbs_bind_surfaces<double,3>(m);

        gbs::ax2<double,3> ax2_z {
                gbs::point<double,3>{0., 0., 0.}, 
                gbs::point<double,3>{0., 0., 1.}, 
                gbs::point<double,3>{1., 0., 0.}
        };
        // Revolution has sense only for 3d case
        py::class_<gbs::SurfaceOfRevolution<double>, std::shared_ptr<gbs::SurfaceOfRevolution<double>>, gbs::Surface<double, 3>>(m, "SurfaceOfRevolution")
        .def(
            py::init<
                std::shared_ptr<gbs::Curve<double, 2>> &,
                const gbs::ax2<double, 3>,
                double, 
                double
            >(),
            py::arg("crv"),
            py::arg("ax") = ax2_z,
            py::arg("a1") = 0.,
            py::arg("a2") = 2 * std::numbers::pi
        )
        .def("__copy__", [](const gbs::SurfaceOfRevolution<double> &self)
                 { return gbs::SurfaceOfRevolution<double>(self); })
        .def(
            py::pickle(
                    [](const gbs::SurfaceOfRevolution<double> &srf) {
                        auto [u1, u2, v1, v2] = srf.bounds();
                            return py::make_tuple(// __getstate__
                                    srf.basisCurve(),
                                    srf.axis2(),
                                    v1,
                                    v2
                            );
                    },
                    [](py::tuple t) { // __setstate__
                            if (t.size() != 4)
                                    throw std::runtime_error("Invalid state!");
                            return gbs::SurfaceOfRevolution<double>{
                                    t[0].cast<std::shared_ptr<gbs::Curve<double, 2>>>(),
                                    t[1].cast<gbs::ax2<double, 3>>(),
                                    t[2].cast<double>(),
                                    t[3].cast<double>(),
                            };
                    }
            )
        )
        .def("__repr__", [](const gbs::SurfaceOfRevolution<double> &self) { return build_rep( self ); } )
        ;
}