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

        py::class_<gbs::BSSfunction<double>>(m,"BSSfunction")
                .def(py::init<const gbs::BSSfunction<double> &>())
                .def(py::init<const gbs::BSSurface<double, 1> &>())
                .def(py::init<const std::vector<double> &, const std::vector<double> &, const std::vector<double> &, size_t, size_t>())
                .def(py::init<const std::vector<double> &, const std::vector<double> &, const std::vector<double> &, const std::vector<size_t> &, const std::vector<size_t> &, size_t, size_t>())
                .def("basisSurface",&gbs::BSSfunction<double>::basisSurface)
                .def("bounds",&gbs::BSSfunction<double>::bounds )
                .def("__call__",&gbs::BSSfunction<double>::operator(),"Function evaluation at given parameter",py::arg("u"),py::arg("v"),py::arg("du") = 0,py::arg("dv") = 0)
                .def("__copy__",  [](const  gbs::BSSfunction<double> &self) {
                        return  gbs::BSSfunction<double>(self);
                })
                .def(py::pickle(
                        [](const gbs::BSSfunction<double> &f) {
                                auto srf = f.basisSurface();
                                return py::make_tuple(
                                                srf.poles(),
                                                srf.knotsU(),
                                                srf.knotsV(),
                                                srf.multsU(),
                                                srf.multsV(),
                                                srf.degreeU(),
                                                srf.degreeV()
                                        );
                                },
                        [](py::tuple t){
                                if (t.size() != 7)
                                        throw std::runtime_error("Invalid state!");
                                gbs::BSSurface<double,1> srf{
                                        t[0].cast<std::vector<std::array<double,1>>>(),
                                        t[1].cast<std::vector<double>>(),
                                        t[2].cast<std::vector<double>>(),
                                        t[3].cast<std::vector<size_t>>(),
                                        t[4].cast<std::vector<size_t>>(),
                                        t[5].cast<size_t>(),
                                        t[6].cast<size_t>()
                                };
                                return gbs::BSSfunction<double>{srf};
                        }
                ))
                // .def("__repr__", [](const gbs::BSSfunction<double> &self) { return build_rep( self ); } )
                ;
}