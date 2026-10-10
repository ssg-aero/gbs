#include "gbsBindBuildSurfaces.h"
#include <gbs/bselementary.h>

namespace
{
    using T = double;

    void gbs_bind_build_elementary(py::module &m)
    {
        using namespace gbs;
        const ax2<T, 3> ax_z{{{0., 0., 0.}, {0., 0., 1.}, {1., 0., 0.}}};
        constexpr T two_pi = 2. * std::numbers::pi;

        m.def("build_circle_arc", py::overload_cast<T, T, T, const ax2<T, 3> &>(&build_circle_arc<T>),
              "Exact arc of circle from theta1 to theta2 (radians) in the frame ax = [center, axis, start direction]; parameter equal to the angle at the knots.",
              py::arg("radius"), py::arg("theta1"), py::arg("theta2"), py::arg("ax") = ax_z);
        m.def("build_circle_arc", py::overload_cast<T, T, T, const point<T, 2> &>(&build_circle_arc<T>),
              "Exact 2D arc of circle from theta1 to theta2 (radians) around center.",
              py::arg("radius"), py::arg("theta1"), py::arg("theta2"), py::arg("center"));
        m.def("build_ellipse_arc", py::overload_cast<T, T, T, T, const ax2<T, 3> &>(&build_ellipse_arc<T>),
              "Exact arc of ellipse O + r1 cos t X + r2 sin t Y, t from theta1 to theta2, in the frame ax = [center, axis, major direction].",
              py::arg("radius1"), py::arg("radius2"), py::arg("theta1"), py::arg("theta2"), py::arg("ax") = ax_z);
        m.def("ellipse_arc_parameter", &ellipse_arc_parameter<T>,
              "Parameter of the point at angle theta of an arc built from theta1 to theta2.",
              py::arg("theta1"), py::arg("theta2"), py::arg("theta"));
        m.def("ellipse_arc_angle", &ellipse_arc_angle<T>,
              "Angle of the point at parameter u of an arc built from theta1 to theta2.",
              py::arg("theta1"), py::arg("theta2"), py::arg("u"));

        auto revolution = [&]<bool rational>() {
            m.def("build_revolution",
                  [](const BSCurveGeneral<T, 3, rational> &g, const ax1<T, 3> &axis, T t1, T t2) { return build_revolution(g, axis, t1, t2); },
                  "Exact rational surface of revolution of a generatrix around axis = [origin, direction]; u is the rotation, v the generatrix parameter.",
                  py::arg("generatrix"), py::arg("axis"), py::arg("theta1") = 0., py::arg("theta2") = two_pi);
        };
        revolution.template operator()<false>();
        revolution.template operator()<true>();
        auto extrusion = [&]<size_t dim, bool rational>() {
            m.def("build_extrusion",
                  [](const BSCurveGeneral<T, dim, rational> &c, const point<T, dim> &V, T v1, T v2) { return build_extrusion(c, V, v1, v2); },
                  "Exact linear extrusion S(u, v) = C(u) + v V, v from v1 to v2.",
                  py::arg("crv"), py::arg("V"), py::arg("v1") = 0., py::arg("v2") = 1.);
        };
        extrusion.template operator()<2, false>();
        extrusion.template operator()<2, true>();
        extrusion.template operator()<3, false>();
        extrusion.template operator()<3, true>();

        m.def("build_cylinder", &build_cylinder<T>, "Exact cylinder O + R (cos u X + sin u Y) + v Z.",
              py::arg("R"), py::arg("ax"), py::arg("v1"), py::arg("v2"), py::arg("theta1") = 0., py::arg("theta2") = two_pi);
        m.def("build_cone", &build_cone<T>, "Exact cone O + (R + v tan(semi_angle)) (cos u X + sin u Y) + v Z.",
              py::arg("R"), py::arg("semi_angle"), py::arg("ax"), py::arg("v1"), py::arg("v2"), py::arg("theta1") = 0., py::arg("theta2") = two_pi);
        m.def("build_sphere", &build_sphere<T>, "Exact sphere O + R cos v (cos u X + sin u Y) + R sin v Z.",
              py::arg("R"), py::arg("ax") = ax_z, py::arg("theta1") = 0., py::arg("theta2") = two_pi,
              py::arg("v1") = -0.5 * std::numbers::pi, py::arg("v2") = 0.5 * std::numbers::pi);
        m.def("build_torus", &build_torus<T>, "Exact torus O + (R + r cos v) (cos u X + sin u Y) + r sin v Z.",
              py::arg("R"), py::arg("r"), py::arg("ax") = ax_z, py::arg("theta1") = 0., py::arg("theta2") = two_pi,
              py::arg("v1") = 0., py::arg("v2") = two_pi);
    }
}

void gbs_bind_build_surface(py::module &m)
{
        gbs_bind_build_surface<double,1>(m);
        gbs_bind_build_surface<double,2>(m);
        gbs_bind_build_surface<double,3>(m);
        gbs_bind_build_elementary(m);
}
