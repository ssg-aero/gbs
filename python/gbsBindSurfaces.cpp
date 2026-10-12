#include "gbsBindSurfaces.h"
#include <gbs/primitives.h>

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
                .def(py::init<const std::vector<double> &, const std::vector<double> &, const std::vector<double> &, size_t, size_t>(),
                        py::arg("poles"), py::arg("knotsFlatsU") , py::arg("knotsFlatsV"), py::arg("degreeU"), py::arg("degreeV")
                )
                .def(py::init<const std::vector<double> &, const std::vector<double> &, const std::vector<double> &, const std::vector<size_t> &, const std::vector<size_t> &, size_t, size_t>(),
                       py::arg("poles"), py::arg("knotsU") , py::arg("knotsV"), py::arg("multsU"), py::arg("multsV"), py::arg("degreeU"), py::arg("degreeV") 
                )
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

        // ---- analytic primitives (gbs/primitives.h), STEP parametrizations on a finite domain
        {
                using namespace gbs;
                using T = double;
                using Box = std::array<T, 4>;
                const Box unit_box{0., 1., 0., 1.};
                const T two_pi = 2. * std::numbers::pi;
                py::class_<ElementaryCurve<T>, std::shared_ptr<ElementaryCurve<T>>, Curve<T, 3>>(m, "ElementaryCurve")
                        .def("parameter", &ElementaryCurve<T>::parameter, py::arg("point"), "Parameter of the nearest point (closed form)")
                        .def("to_nurbs", &ElementaryCurve<T>::to_nurbs, "Exact NURBS over the range")
                        .def("nurbs_parameter", &ElementaryCurve<T>::nurbs_parameter, py::arg("t"))
                        .def("reversed", &ElementaryCurve<T>::reversed, "Same points run backwards: t -> -t");
                py::class_<LineCurve<T>, std::shared_ptr<LineCurve<T>>, ElementaryCurve<T>>(m, "LineCurve", "C(t) = O + t V, t in [t1, t2]")
                        .def(py::init<const point<T, 3> &, const point<T, 3> &, T, T>(), py::arg("origin"), py::arg("vector"), py::arg("t1"), py::arg("t2"))
                        .def("origin", &LineCurve<T>::origin)
                        .def("vector", &LineCurve<T>::vector);
                py::class_<EllipseCurve<T>, std::shared_ptr<EllipseCurve<T>>, ElementaryCurve<T>>(m, "EllipseCurve",
                                                                                               "C(t) = O + r1 cos t X + r2 sin t Y in the frame ax = [center, axis, major direction]")
                        .def(py::init<T, T, const ax2<T, 3> &, T, T>(), py::arg("r1"), py::arg("r2"), py::arg("ax"), py::arg("t1") = 0., py::arg("t2") = two_pi)
                        .def("radius1", &EllipseCurve<T>::radius1)
                        .def("radius2", &EllipseCurve<T>::radius2)
                        .def("is_circle", &EllipseCurve<T>::is_circle);

                py::class_<ElementarySurface<T>, std::shared_ptr<ElementarySurface<T>>, Surface<T, 3>>(m, "ElementarySurface")
                        .def("parameters", &ElementarySurface<T>::parameters, py::arg("point"), "Parameters of a point (closed form)")
                        .def("to_nurbs", &ElementarySurface<T>::to_nurbs, "Exact NURBS over the box")
                        .def("nurbs_parameters", &ElementarySurface<T>::nurbs_parameters, py::arg("u"), py::arg("v"))
                        .def("iso_u", &ElementarySurface<T>::iso_u, py::arg("u"), py::arg("v1"), py::arg("v2"));
                py::class_<Plane<T>, std::shared_ptr<Plane<T>>, ElementarySurface<T>>(m, "Plane", "S(u, v) = O + u X + v Y")
                        .def(py::init<const ax2<T, 3> &, Box>(), py::arg("ax") = ax2_z, py::arg("box") = unit_box);
                py::class_<CylindricalSurface<T>, std::shared_ptr<CylindricalSurface<T>>, ElementarySurface<T>>(m, "CylindricalSurface",
                                                                                                             "S(u, v) = O + R (cos u X + sin u Y) + v Z")
                        .def(py::init<T, const ax2<T, 3> &, Box>(), py::arg("radius"), py::arg("ax") = ax2_z, py::arg("box") = Box{0., two_pi, 0., 1.})
                        .def("radius", &CylindricalSurface<T>::radius);
                py::class_<ConicalSurface<T>, std::shared_ptr<ConicalSurface<T>>, ElementarySurface<T>>(m, "ConicalSurface",
                                                                                                     "S(u, v) = O + (R + v tan a)(cos u X + sin u Y) + v Z")
                        .def(py::init<T, T, const ax2<T, 3> &, Box>(), py::arg("radius"), py::arg("semi_angle"), py::arg("ax") = ax2_z, py::arg("box") = Box{0., two_pi, 0., 1.})
                        .def("radius", &ConicalSurface<T>::radius)
                        .def("semi_angle", &ConicalSurface<T>::semi_angle);
                py::class_<SphericalSurface<T>, std::shared_ptr<SphericalSurface<T>>, ElementarySurface<T>>(m, "SphericalSurface",
                                                                                                         "S(u, v) = O + R cos v (cos u X + sin u Y) + R sin v Z")
                        .def(py::init<T, const ax2<T, 3> &, Box>(), py::arg("radius"), py::arg("ax") = ax2_z,
                             py::arg("box") = Box{0., two_pi, -0.5 * std::numbers::pi, 0.5 * std::numbers::pi})
                        .def("radius", &SphericalSurface<T>::radius);
                py::class_<ToroidalSurface<T>, std::shared_ptr<ToroidalSurface<T>>, ElementarySurface<T>>(m, "ToroidalSurface",
                                                                                                       "S(u, v) = O + (R + r cos v)(cos u X + sin u Y) + r sin v Z")
                        .def(py::init<T, T, const ax2<T, 3> &, Box>(), py::arg("major"), py::arg("minor"), py::arg("ax") = ax2_z, py::arg("box") = Box{0., two_pi, 0., two_pi})
                        .def("major_radius", &ToroidalSurface<T>::major_radius)
                        .def("minor_radius", &ToroidalSurface<T>::minor_radius);
        }
}
