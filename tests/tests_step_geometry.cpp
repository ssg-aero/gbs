#include <doctest_gtest.hpp>
#include <gbs-io/step/geometry.h>

#include <cmath>
#include <functional>
#include <numbers>
#include <string>

#ifdef GBS_USE_MODULES
    import vecop;
#endif

using namespace gbs;
using namespace gbs::step;
using gbs::operator-;
using gbs::operator+;
using gbs::operator*;

namespace
{
    constexpr double pi = std::numbers::pi;
    constexpr double tol = 1e-12;

    // A STEP file around DATA lines, with a context in millimetres and radians as #1.
    std::string step_file(const std::string &data, const std::string &units = "")
    {
        const std::string ctx = units.empty() ? R"(
#1=(GEOMETRIC_REPRESENTATION_CONTEXT(3)GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT((#5))GLOBAL_UNIT_ASSIGNED_CONTEXT((#2,#3,#4))REPRESENTATION_CONTEXT('',''));
#2=(LENGTH_UNIT()NAMED_UNIT(*)SI_UNIT(.MILLI.,.METRE.));
#3=(NAMED_UNIT(*)PLANE_ANGLE_UNIT()SI_UNIT($,.RADIAN.));
#4=(NAMED_UNIT(*)SI_UNIT($,.STERADIAN.)SOLID_ANGLE_UNIT());
#5=UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(1.E-07),#2,'distance_accuracy_value','');
)" : units;
        return "ISO-10303-21;\nHEADER;\nFILE_DESCRIPTION((''),'2;1');\nFILE_NAME('','',(''),(''),'','','');\n"
               "FILE_SCHEMA(('AUTOMOTIVE_DESIGN'));\nENDSEC;\nDATA;\n" +
               ctx + data + "ENDSEC;\nEND-ISO-10303-21;\n";
    }

    P21File parse(const std::string &s)
    {
        auto f = parse_p21(s);
        if (!f)
            throw std::runtime_error(f.error().message);
        return std::move(*f);
    }

    double dist(const Point<3> &a, const Point<3> &b) { return norm(a - b); }

    // Frame used by most entities: origin (1,2,3), axis (0,0,1), reference (1,0,0)
    const std::string frame = R"(
#10=CARTESIAN_POINT('',(1.,2.,3.));
#11=DIRECTION('',(0.,0.,1.));
#12=DIRECTION('',(1.,0.,0.));
#13=AXIS2_PLACEMENT_3D('',#10,#11,#12);
)";

    // Cox-de Boor evaluation of a B-spline on any knot vector, for reference.
    Point<3> de_boor(const std::vector<double> &k, const std::vector<Point<3>> &P, std::size_t p, double u)
    {
        std::size_t s = p;
        while (s + 1 < k.size() - p - 1 && u >= k[s + 1])
            ++s;
        std::vector<Point<3>> d(P.begin() + static_cast<long>(s - p), P.begin() + static_cast<long>(s + 1));
        for (std::size_t r = 1; r <= p; ++r)
            for (std::size_t j = p; j >= r; --j)
            {
                const std::size_t i = j + s - p;
                const double a = (u - k[i]) / (k[i + p + 1 - r] - k[i]);
                d[j] = (1. - a) * d[j - 1] + a * d[j];
            }
        return d[p];
    }
}

TEST(tests_step_geometry, units)
{
    // millimetres, radians, uncertainty
    {
        auto f = parse(step_file(""));
        auto u = read_units(f, 1);
        ASSERT_NEAR(u.length_mm, 1., 1e-15);
        ASSERT_EQ(u.length, 1.);
        ASSERT_EQ(u.angle, 1.);
        ASSERT_EQ(u.length_name, "MILLI METRE");
        ASSERT_EQ(u.angle_name, "RADIAN");
        ASSERT_NEAR(*u.uncertainty, 1e-7, 1e-20);
        auto in_inch = read_units(f, 1, 25.4);
        ASSERT_NEAR(in_inch.length, 1. / 25.4, 1e-15);
        ASSERT_NEAR(*in_inch.uncertainty, 1e-7 / 25.4, 1e-20);
    }
    // inches and degrees as conversion based units, simple instances
    {
        auto f = parse(step_file("", R"(
#1=(GEOMETRIC_REPRESENTATION_CONTEXT(3)GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT((#9))GLOBAL_UNIT_ASSIGNED_CONTEXT((#2,#5))REPRESENTATION_CONTEXT('',''));
#2=(CONVERSION_BASED_UNIT('INCH',#3)LENGTH_UNIT()NAMED_UNIT(#4));
#3=LENGTH_MEASURE_WITH_UNIT(LENGTH_MEASURE(25.4),#6);
#4=DIMENSIONAL_EXPONENTS(1.,0.,0.,0.,0.,0.,0.);
#5=(CONVERSION_BASED_UNIT('DEGREE',#7)NAMED_UNIT(#4)PLANE_ANGLE_UNIT());
#6=(LENGTH_UNIT()NAMED_UNIT(*)SI_UNIT(.MILLI.,.METRE.));
#7=(MEASURE_WITH_UNIT(PLANE_ANGLE_MEASURE(0.0174532925199433),#8)PLANE_ANGLE_MEASURE_WITH_UNIT());
#8=(NAMED_UNIT(*)PLANE_ANGLE_UNIT()SI_UNIT($,.RADIAN.));
#9=UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(0.001),#2,'distance_accuracy_value','');
)"));
        auto u = read_units(f, 1);
        ASSERT_NEAR(u.length_mm, 25.4, 1e-12);
        ASSERT_NEAR(u.length, 25.4, 1e-12);
        ASSERT_NEAR(u.angle, pi / 180., 1e-15);
        ASSERT_EQ(u.length_name, "INCH");
        ASSERT_EQ(u.angle_name, "DEGREE");
        ASSERT_NEAR(*u.uncertainty, 0.0254, 1e-15);
    }
    // metres as a simple SI_UNIT, no uncertainty
    {
        auto f = parse(step_file("", R"(
#1=(GEOMETRIC_REPRESENTATION_CONTEXT(3)GLOBAL_UNIT_ASSIGNED_CONTEXT((#2))REPRESENTATION_CONTEXT('',''));
#2=(LENGTH_UNIT()NAMED_UNIT(*)SI_UNIT($,.METRE.));
)"));
        auto u = read_units(f, 1);
        ASSERT_EQ(u.length, 1000.);
        ASSERT_FALSE(u.uncertainty.has_value());
    }
    // unknown unit
    {
        auto f = parse(step_file("", R"(
#1=(GEOMETRIC_REPRESENTATION_CONTEXT(3)GLOBAL_UNIT_ASSIGNED_CONTEXT((#2))REPRESENTATION_CONTEXT('',''));
#2=(LENGTH_UNIT()NAMED_UNIT(*)SI_UNIT($,.FOOT.));
)"));
        ASSERT_THROW(read_units(f, 1), StepError);
    }
}

TEST(tests_step_geometry, lines_and_conics)
{
    auto f = parse(step_file(frame + R"(
#20=DIRECTION('',(0.,1.,0.));
#21=VECTOR('',#20,2.);
#22=LINE('',#10,#21);
#30=CIRCLE('',#13,5.);
#31=ELLIPSE('',#13,4.,1.5);
#40=CARTESIAN_POINT('',(0.,0.));
#41=DIRECTION('',(0.,1.));
#42=AXIS2_PLACEMENT_2D('',#40,#41);
#43=CIRCLE('',#42,0.5);
)"));
    // in millimetres, then with lengths multiplied by 1000 (a file in metres)
    for (double L : {1., 1000.})
    {
        GeometryReader r{f, Units{.length = L}};
        const Point<3> O{L * 1., L * 2., L * 3.};
        // line: C(t) = O + t (0, 2, 0), t unchanged by the unit
        auto line = r.curve<3>(22);
        ASSERT_LT(dist(line->value(1.5), O + Point<3>{0., 3. * L, 0.}), tol * L);
        ASSERT_NEAR(line->parameter(line->value(-0.7)), -0.7, tol);
        ASSERT_FALSE(line->periodic());
        auto seg = line->to_nurbs(-1., 2.);
        ASSERT_EQ(seg.u1, -1.);
        ASSERT_LT(dist(seg.curve->value(0.5), line->value(0.5)), tol * L);

        // circle and ellipse: STEP parameter is the angle
        auto c = r.curve<3>(30);
        auto e = r.curve<3>(31);
        ASSERT_TRUE(c->periodic());
        ASSERT_EQ(c->kind(), ParamKind::Angle);
        for (double t : {0., 1., 2.5, 4., 6.})
        {
            ASSERT_LT(dist(c->value(t), O + Point<3>{5. * L * std::cos(t), 5. * L * std::sin(t), 0.}), 1e-12 * L);
            ASSERT_LT(dist(e->value(t), O + Point<3>{4. * L * std::cos(t), 1.5 * L * std::sin(t), 0.}), 1e-12 * L);
            ASSERT_NEAR(c->parameter(c->value(t)), t, 1e-12);
            ASSERT_NEAR(e->parameter(e->value(t)), t, 1e-12);
        }
        // arcs from 5 to 5 + 3 (across 2 pi): exact NURBS, ends at the requested parameters
        for (auto *d : {c.get(), e.get()})
        {
            auto n = d->to_nurbs(5., 8.);
            ASSERT_EQ(n.u1, 5.);
            ASSERT_EQ(n.u2, 8.);
            for (double t : {5., 6., 7.3, 8.})
                ASSERT_LT(dist(n.curve->value(d->nurbs_parameter(5., 8., t)), d->value(t)), 1e-11 * L);
        }
        ASSERT_THROW(c->to_nurbs(0., 7.), StepError);
    }
    // 2D circle: parameter space, no unit conversion
    GeometryReader r{f, Units{.length = 1000.}};
    auto c2 = r.curve<2>(43);
    ASSERT_NEAR(c2->value(0.)[1], 0.5, tol); // reference direction (0, 1)
    ASSERT_NEAR(c2->value(pi / 2)[0], -0.5, tol);
}

TEST(tests_step_geometry, bsplines)
{
    // the same quadratic curve written four ways; a rational quarter circle; a periodic (unclamped) curve
    auto f = parse(step_file(R"(
#10=CARTESIAN_POINT('',(0.,0.,0.));
#11=CARTESIAN_POINT('',(1.,2.,0.));
#12=CARTESIAN_POINT('',(3.,2.,1.));
#13=CARTESIAN_POINT('',(4.,0.,0.));
#14=CARTESIAN_POINT('',(5.,-1.,2.));
#20=B_SPLINE_CURVE_WITH_KNOTS('',2,(#10,#11,#12,#13,#14),.UNSPECIFIED.,.F.,.F.,(3,1,1,3),(0.,1.,2.,3.),.UNSPECIFIED.);
#21=QUASI_UNIFORM_CURVE('',2,(#10,#11,#12,#13,#14),.UNSPECIFIED.,.F.,.F.);
#22=BEZIER_CURVE('',2,(#10,#11,#12,#13,#14),.UNSPECIFIED.,.F.,.F.);
#23=UNIFORM_CURVE('',2,(#10,#11,#12,#13,#14),.UNSPECIFIED.,.F.,.F.);
#30=CARTESIAN_POINT('',(1.,0.,0.));
#31=CARTESIAN_POINT('',(1.,1.,0.));
#32=CARTESIAN_POINT('',(0.,1.,0.));
#33=(BOUNDED_CURVE()B_SPLINE_CURVE(2,(#30,#31,#32),.CIRCULAR_ARC.,.F.,.F.)B_SPLINE_CURVE_WITH_KNOTS((3,3),(0.,1.),.UNSPECIFIED.)CURVE()GEOMETRIC_REPRESENTATION_ITEM()RATIONAL_B_SPLINE_CURVE((1.,0.707106781186548,1.))REPRESENTATION_ITEM(''));
#40=B_SPLINE_CURVE_WITH_KNOTS('',3,(#10,#11,#12,#13,#14,#10,#11,#12),.UNSPECIFIED.,.T.,.F.,(1,1,1,1,1,1,1,1,1,1,1,1),(-3.,-2.,-1.,0.,1.,2.,3.,4.,5.,6.,7.,8.),.UNSPECIFIED.);
)"));
    GeometryReader r{f, Units{.length = 2.}};
    const std::vector<Point<3>> P{{0., 0., 0.}, {2., 4., 0.}, {6., 4., 2.}, {8., 0., 0.}, {10., -2., 4.}}; // scaled by 2
    // with knots and quasi uniform: same knots (0 0 0 1 2 3 3 3)
    for (std::uint64_t id : {20, 21})
    {
        auto c = r.curve<3>(id);
        ASSERT_EQ(c->domain()[0], 0.);
        ASSERT_EQ(c->domain()[1], 3.);
        for (double u : {0., 0.4, 1., 1.7, 2.9, 3.})
            ASSERT_LT(dist(c->value(u), de_boor({0, 0, 0, 1, 2, 3, 3, 3}, P, 2, u)), tol);
        ASSERT_NEAR(c->parameter(c->value(1.3)), 1.3, 1e-10);
    }
    // Bezier: two quadratic segments on (0 0 0 1 1 2 2 2)
    {
        auto c = r.curve<3>(22);
        for (double u : {0., 0.5, 1., 1.5, 2.})
            ASSERT_LT(dist(c->value(u), de_boor({0, 0, 0, 1, 1, 2, 2, 2}, P, 2, u)), tol);
    }
    // uniform: unclamped knots (-2 … 5), domain [0, 3], clamped by knot insertion
    {
        auto c = r.curve<3>(23);
        ASSERT_NEAR(c->domain()[0], 0., tol);
        ASSERT_NEAR(c->domain()[1], 3., tol);
        for (double u : {0., 0.3, 1., 2.2, 3.})
            ASSERT_LT(dist(c->value(u), de_boor({-2, -1, 0, 1, 2, 3, 4, 5}, P, 2, u)), 1e-12);
    }
    // rational quarter circle of radius 2 (scaled): homogeneous poles built from Cartesian points and weights
    {
        auto c = r.curve<3>(33);
        for (int i = 0; i <= 10; ++i)
            ASSERT_NEAR(norm(c->value(i / 10.)), 2., 1e-12);
        ASSERT_LT(dist(c->value(0.5), Point<3>{std::sqrt(2.), std::sqrt(2.), 0.}), 1e-12);
    }
    // periodic cubic written with unclamped knots: domain [0, 5]
    {
        const std::vector<double> k{-3., -2., -1., 0., 1., 2., 3., 4., 5., 6., 7., 8.};
        std::vector<Point<3>> Q{P[0], P[1], P[2], P[3], P[4], P[0], P[1], P[2]};
        auto c = r.curve<3>(40);
        ASSERT_NEAR(c->domain()[0], 0., tol);
        ASSERT_NEAR(c->domain()[1], 5., tol);
        for (double u : {0., 0.7, 2.5, 4.1, 5.})
            ASSERT_LT(dist(c->value(u), de_boor(k, Q, 3, u)), 1e-12);
        ASSERT_LT(dist(c->value(0.), c->value(5.)), 1e-12); // closed
    }
}

TEST(tests_step_geometry, trimmed_composite_and_polyline)
{
    auto f = parse(step_file(frame + R"(
#30=CIRCLE('',#13,2.);
#31=CARTESIAN_POINT('',(1.,4.,3.));
#40=TRIMMED_CURVE('',#30,(PARAMETER_VALUE(300.)),(PARAMETER_VALUE(60.)),.T.,.PARAMETER.);
#41=TRIMMED_CURVE('',#30,(PARAMETER_VALUE(300.)),(PARAMETER_VALUE(60.)),.F.,.PARAMETER.);
#42=TRIMMED_CURVE('',#30,(#31,PARAMETER_VALUE(10.)),(PARAMETER_VALUE(180.)),.T.,.CARTESIAN.);
#50=CARTESIAN_POINT('',(0.,0.,0.));
#51=CARTESIAN_POINT('',(1.,0.,0.));
#52=CARTESIAN_POINT('',(1.,1.,0.));
#53=POLYLINE('',(#50,#51,#52));
#60=DIRECTION('',(1.,0.,0.));
#61=VECTOR('',#60,1.);
#62=LINE('',#50,#61);
#63=TRIMMED_CURVE('',#62,(PARAMETER_VALUE(0.)),(PARAMETER_VALUE(1.)),.T.,.PARAMETER.);
#64=DIRECTION('',(0.,1.,0.));
#65=VECTOR('',#64,1.);
#66=LINE('',#51,#65);
#67=TRIMMED_CURVE('',#66,(PARAMETER_VALUE(1.)),(PARAMETER_VALUE(0.)),.F.,.PARAMETER.);
#68=COMPOSITE_CURVE_SEGMENT(.CONTINUOUS.,.T.,#63);
#69=COMPOSITE_CURVE_SEGMENT(.CONTINUOUS.,.F.,#67);
#70=COMPOSITE_CURVE('',(#68,#69),.F.);
#80=SURFACE_CURVE('',#30,(),.CURVE_3D.);
#81=HYPERBOLA('',#13,1.,1.);
)"));
    GeometryReader r{f, Units{.angle = pi / 180.}};
    // degrees: from 300 to 60 forward crosses 0 -> [300, 420] degrees
    auto t1 = r.curve<3>(40);
    ASSERT_NEAR(t1->domain()[0], 300. * pi / 180., 1e-12);
    ASSERT_NEAR(t1->domain()[1], 420. * pi / 180., 1e-12);
    ASSERT_FALSE(t1->reversed());
    ASSERT_NEAR(t1->parameter(t1->value(0.2)), 2 * pi + 0.2, 1e-12); // the turn inside the trim
    // backward from 300 to 60: the other side, [60, 300], reversed
    auto t2 = r.curve<3>(41);
    ASSERT_NEAR(t2->domain()[0], pi / 3, 1e-12);
    ASSERT_NEAR(t2->domain()[1], 5 * pi / 3, 1e-12);
    ASSERT_TRUE(t2->reversed());
    // trimmed by a point (master representation CARTESIAN): (1, 4, 3) is at 90 degrees
    auto t3 = r.curve<3>(42);
    ASSERT_NEAR(t3->domain()[0], pi / 2, 1e-12);
    ASSERT_NEAR(t3->domain()[1], pi, 1e-12);
    auto n3 = t3->to_nurbs(pi / 2, pi);
    ASSERT_LT(dist(n3.curve->begin(), Point<3>{1., 4., 3.}), 1e-12);
    ASSERT_LT(dist(n3.curve->end(), Point<3>{-1., 2., 3.}), 1e-12);

    // polyline: one unit of parameter per segment
    auto pl = r.curve<3>(53);
    ASSERT_EQ(pl->domain()[1], 2.);
    ASSERT_LT(dist(pl->value(1.5), Point<3>{1., 0.5, 0.}), tol);
    // composite: segment (0,0,0)->(1,0,0), then the reversed of a reversed trim (1,0,0)->(1,1,0)
    auto cc = r.curve<3>(70);
    const auto [a, b] = cc->domain();
    ASSERT_LT(dist(cc->value(a), Point<3>{0., 0., 0.}), tol);
    ASSERT_LT(dist(cc->value(b), Point<3>{1., 1., 0.}), tol);
    ASSERT_LT(dist(cc->value(0.5 * (a + b)), Point<3>{1., 0., 0.}), 1e-12);
    // surface curve: its 3D curve
    ASSERT_LT(dist(r.curve<3>(80)->value(0.), Point<3>{3., 2., 3.}), tol);
    // unsupported
    try
    {
        r.curve<3>(81);
        FAIL("HYPERBOLA should be unsupported");
    }
    catch (const StepUnsupported &e)
    {
        ASSERT_EQ(e.type(), "HYPERBOLA");
        ASSERT_EQ(e.id(), 81u);
    }
}

TEST(tests_step_geometry, elementary_surfaces)
{
    auto f = parse(step_file(frame + R"(
#20=PLANE('',#13);
#21=CYLINDRICAL_SURFACE('',#13,2.);
#22=CONICAL_SURFACE('',#13,1.,30.);
#23=SPHERICAL_SURFACE('',#13,3.);
#24=TOROIDAL_SURFACE('',#13,4.,1.);
#25=DEGENERATE_TOROIDAL_SURFACE('',#13,0.5,1.,.T.);
#26=OFFSET_SURFACE('',#21,0.5,.F.);
#27=OFFSET_SURFACE('',#22,0.2,.F.);
#28=RECTANGULAR_TRIMMED_SURFACE('',#21,-30.,90.,0.,5.,.T.,.T.);
)"));
    const double L = 10.; // centimetres
    GeometryReader r{f, Units{.length = L, .angle = pi / 180.}};
    const Point<3> O{L, 2. * L, 3. * L};
    const Point<3> X{1., 0., 0.}, Y{0., 1., 0.}, Z{0., 0., 1.};
    auto rad = [&](double u) { return std::cos(u) * X + std::sin(u) * Y; };
    const double ta = std::tan(pi / 6);

    struct Case
    {
        std::uint64_t id;
        std::function<Point<3>(double, double)> exact;
        UVBox box;
    };
    const std::vector<Case> cases{
        {20, [&](double u, double v) { return O + u * X + v * Y; }, {-5., 7., -3., 2.}},
        {21, [&](double u, double v) { return O + 2. * L * rad(u) + v * Z; }, {-0.5, 2., -4., 9.}},
        {22, [&](double u, double v) { return O + (L + v * ta) * rad(u) + v * Z; }, {0., 2 * pi, -L / ta, 20.}},
        {23, [&](double u, double v) { return O + 3. * L * std::cos(v) * rad(u) + 3. * L * std::sin(v) * Z; }, {0., 2 * pi, -pi / 2, pi / 2}},
        {24, [&](double u, double v) { return O + (4. * L + L * std::cos(v)) * rad(u) + L * std::sin(v) * Z; }, {1., 3., 0., 2 * pi}},
        {25, [&](double u, double v) { return O + (0.5 * L + L * std::cos(v)) * rad(u) + L * std::sin(v) * Z; }, {0., 2 * pi, -1., 1.}},
        {26, [&](double u, double v) { return O + 2.5 * L * rad(u) + v * Z; }, {0., 2 * pi, 0., 5.}},
        {28, [&](double u, double v) { return O + 2. * L * rad(u) + v * Z; }, {-pi / 6, pi / 2, 0., 5. * L}},
    };
    for (const auto &c : cases)
    {
        auto s = r.surface(c.id);
        auto nurbs = s->to_nurbs(c.box);
        for (int i = 0; i <= 8; ++i)
            for (int j = 0; j <= 8; ++j)
            {
                const double u = c.box.u1 + (c.box.u2 - c.box.u1) * i / 8., v = c.box.v1 + (c.box.v2 - c.box.v1) * j / 8.;
                const auto p = c.exact(u, v);
                ASSERT_LT(dist(s->value(u, v), p), 1e-12 * L);
                auto [nu, nv] = s->nurbs_parameters(c.box, u, v);
                ASSERT_LT(dist(nurbs->value(nu, nv), p), 1e-11 * L);
                if (std::abs(std::cos(v)) > 1e-6 || c.id != 23) // a sphere pole has no u
                {
                    auto uv = *s->parameters(p);
                    ASSERT_LT(dist(s->value(uv[0], uv[1]), p), 1e-11 * L);
                }
            }
    }
    // the trimmed surface keeps its trims, converted: degrees for u, lengths for v
    auto tr = r.surface(28)->domain();
    ASSERT_NEAR(tr.u1, -pi / 6, 1e-15);
    ASSERT_NEAR(tr.v2, 5. * L, 1e-12);
    // offset cone: points at distance 0.2 L from the cone, along its normal
    auto cone = r.surface(22);
    auto off = r.surface(27);
    for (double v : {-5., 0., 7.})
    {
        const double u = 0.8;
        const auto n = std::cos(pi / 6) * rad(u) - std::sin(pi / 6) * Z;
        const auto p = cone->value(u, v) + 0.2 * L * n;
        auto uv = *off->parameters(p);
        ASSERT_LT(dist(off->value(uv[0], uv[1]), p), 1e-11 * L);
    }
}

TEST(tests_step_geometry, swept_and_bspline_surfaces)
{
    auto f = parse(step_file(frame + R"(
#30=CIRCLE('',#13,1.);
#40=CARTESIAN_POINT('',(4.,2.,3.));
#41=DIRECTION('',(1.,0.,1.));
#42=VECTOR('',#41,1.);
#43=LINE('',#40,#42);
#44=AXIS1_PLACEMENT('',#10,#11);
#45=SURFACE_OF_REVOLUTION('',#43,#44);
#46=DIRECTION('',(0.,1.,0.));
#47=AXIS2_PLACEMENT_3D('',#40,#46,#12);
#48=CIRCLE('',#47,0.5);
#49=SURFACE_OF_REVOLUTION('',#48,#44);
#50=DIRECTION('',(0.,1.,1.));
#51=VECTOR('',#50,3.);
#52=SURFACE_OF_LINEAR_EXTRUSION('',#30,#51);
#53=SURFACE_OF_LINEAR_EXTRUSION('',#43,#51);
#60=CARTESIAN_POINT('',(0.,0.,0.));
#61=CARTESIAN_POINT('',(0.,1.,0.));
#62=CARTESIAN_POINT('',(0.,2.,1.));
#63=CARTESIAN_POINT('',(1.,0.,0.));
#64=CARTESIAN_POINT('',(1.,1.,2.));
#65=CARTESIAN_POINT('',(1.,2.,0.));
#66=B_SPLINE_SURFACE_WITH_KNOTS('',1,2,((#60,#61,#62),(#63,#64,#65)),.UNSPECIFIED.,.F.,.F.,.F.,(2,2),(3,3),(0.,1.),(0.,2.),.UNSPECIFIED.);
#67=(BOUNDED_SURFACE()B_SPLINE_SURFACE(1,2,((#60,#61,#62),(#63,#64,#65)),.UNSPECIFIED.,.F.,.F.,.F.)B_SPLINE_SURFACE_WITH_KNOTS((2,2),(3,3),(0.,1.),(0.,2.),.UNSPECIFIED.)GEOMETRIC_REPRESENTATION_ITEM()RATIONAL_B_SPLINE_SURFACE(((1.,2.,1.),(1.,2.,1.)))REPRESENTATION_ITEM('')SURFACE());
)"));
    GeometryReader r{f, Units{}};
    const Point<3> O{1., 2., 3.}, Z{0., 0., 1.};
    auto rotate = [&](const Point<3> &p, double a) {
        const auto q = p - O;
        const auto h = (q * Z) * Z;
        return O + h + std::cos(a) * (q - h) + std::sin(a) * cross(Z, q - h);
    };
    // revolution of a line (a cone) and of a circle (a torus around a tilted axis frame)
    for (auto [id, gid, box] : std::vector<std::tuple<std::uint64_t, std::uint64_t, UVBox>>{
             {45, 43, {0., 2 * pi, -1., 2.}}, {49, 48, {0.5, 2., 0., 2 * pi}}})
    {
        auto s = r.surface(id);
        auto g = r.curve<3>(gid);
        auto n = s->to_nurbs(box);
        for (int i = 0; i <= 6; ++i)
            for (int j = 0; j <= 6; ++j)
            {
                const double u = box.u1 + (box.u2 - box.u1) * i / 6., v = box.v1 + (box.v2 - box.v1) * j / 6.;
                const auto p = rotate(g->value(v), u);
                ASSERT_LT(dist(s->value(u, v), p), 1e-12);
                auto [nu, nv] = s->nurbs_parameters(box, u, v);
                ASSERT_LT(dist(n->value(nu, nv), p), 1e-11);
                auto uv = *s->parameters(p);
                ASSERT_LT(dist(s->value(uv[0], uv[1]), p), 1e-10);
            }
    }
    // extrusions of a circle and of a line along V = 3 (0, 1, 1) / sqrt 2
    const Point<3> V = (3. / std::sqrt(2.)) * Point<3>{0., 1., 1.};
    for (auto [id, cid, box] : std::vector<std::tuple<std::uint64_t, std::uint64_t, UVBox>>{
             {52, 30, {-1., 2., -0.5, 1.5}}, {53, 43, {-2., 1., 0., 2.}}})
    {
        auto s = r.surface(id);
        auto c = r.curve<3>(cid);
        auto n = s->to_nurbs(box);
        for (int i = 0; i <= 6; ++i)
            for (int j = 0; j <= 6; ++j)
            {
                const double u = box.u1 + (box.u2 - box.u1) * i / 6., v = box.v1 + (box.v2 - box.v1) * j / 6.;
                const auto p = c->value(u) + v * V;
                ASSERT_LT(dist(s->value(u, v), p), 1e-12);
                auto [nu, nv] = s->nurbs_parameters(box, u, v);
                ASSERT_LT(dist(n->value(nu, nv), p), 1e-11);
                auto uv = *s->parameters(p);
                ASSERT_LT(dist(s->value(uv[0], uv[1]), p), 1e-10);
                ASSERT_NEAR(uv[1], v, 1e-10);
            }
    }
    // B-spline surfaces: STEP [u][v] grid, gbs u fastest; corners and a rational weight
    auto sp = r.surface(66);
    auto b = sp->domain();
    ASSERT_EQ(b.u2, 1.);
    ASSERT_EQ(b.v2, 2.);
    ASSERT_FALSE(sp->parameters({0., 0., 0.}).has_value());
    ASSERT_LT(dist(sp->value(0., 2.), Point<3>{0., 2., 1.}), tol);
    ASSERT_LT(dist(sp->value(1., 0.), Point<3>{1., 0., 0.}), tol);
    // v = 1 is the middle of the quadratic Bezier in v: rows (0, 1, 0.25) at u = 0 and (1, 1, 1) at u = 1
    ASSERT_LT(dist(sp->value(0.5, 1.), Point<3>{0.5, 1., 0.625}), tol);
    auto rs = r.surface(67);
    // weights (1, 2, 1) along v: weighted mean of (0,0,0), (0,1,0), (0,2,1) with Bernstein (1/4, 1/2, 1/4)
    const auto mid = rs->value(0., 1.);
    ASSERT_LT(dist(mid, Point<3>{0., (0.25 * 0. + 0.5 * 2. * 1. + 0.25 * 2.) / (0.25 + 1. + 0.25), (0.25 * 1.) / 1.5}), tol);
}

TEST(tests_step_geometry, parameter_boxes)
{
    auto f = parse(step_file(frame + R"(
#20=PLANE('',#13);
#21=CYLINDRICAL_SURFACE('',#13,2.);
#23=SPHERICAL_SURFACE('',#13,3.);
#22=CONICAL_SURFACE('',#13,1.,0.5);
)"));
    GeometryReader r{f, Units{}};
    auto cyl = r.surface(21);
    auto loop = [&](const SurfaceDef &s, std::vector<std::array<double, 2>> uv) {
        std::vector<Point<3>> l;
        for (auto [u, v] : uv)
            l.push_back(s.value(u, v));
        return l;
    };
    // a patch of cylinder from -20 to 30 degrees: does not reach the seam, box across 0
    {
        const double a = -pi / 9, b = pi / 6;
        std::vector<std::vector<Point<3>>> loops{loop(*cyl, {{a, 0.}, {0., 0.}, {b, 0.}, {b, 4.}, {0., 4.}, {a, 4.}})};
        auto box = parameter_box(*cyl, loops);
        const double m = 0.1 * (b - a);
        ASSERT_NEAR(box.u1, a - m, 1e-12);
        ASSERT_NEAR(box.u2, b + m, 1e-12);
        ASSERT_NEAR(box.v1, -0.4, 1e-12);
        ASSERT_NEAR(box.v2, 4.4, 1e-12);
    }
    // a band around the cylinder (a loop winding once): full turn, the STEP seam
    {
        std::vector<Point<3>> l;
        for (int i = 0; i < 12; ++i)
            l.push_back(cyl->value(i * pi / 6, 1.));
        std::vector<std::vector<Point<3>>> loops{l};
        auto box = parameter_box(*cyl, loops);
        ASSERT_EQ(box.u1, 0.);
        ASSERT_EQ(box.u2, 2 * pi);
    }
    // two loops seen across the 0 / 2 pi cut are joined, not spread over the turn
    {
        std::vector<std::vector<Point<3>>> loops{loop(*cyl, {{6.0, 0.}, {6.2, 0.}, {6.2, 1.}, {6.0, 1.}}),
                                                 loop(*cyl, {{0.1, 0.}, {0.2, 0.}, {0.2, 1.}})};
        auto box = parameter_box(*cyl, loops);
        ASSERT_LT(box.u2 - box.u1, 0.6);
    }
    // a sphere cap with its pole: latitude clamped to pi / 2, pole ignored in u
    {
        auto sph = r.surface(23);
        std::vector<Point<3>> l;
        for (int i = 0; i < 12; ++i)
            l.push_back(sph->value(i * pi / 6, pi / 4));
        l.push_back(sph->value(0., pi / 2));
        std::vector<std::vector<Point<3>>> loops{l};
        auto box = parameter_box(*sph, loops);
        ASSERT_EQ(box.u2, 2 * pi);
        ASSERT_NEAR(box.v2, pi / 2, 1e-12);
        ASSERT_NEAR(box.v1, pi / 4 - 0.1 * pi / 4, 1e-12);
    }
    // a cone face reaching its apex: v clamped at the apex
    {
        auto cone = r.surface(22);
        const double apex = -1. / std::tan(0.5);
        std::vector<std::vector<Point<3>>> loops{loop(*cone, {{0., apex}, {1., 0.}, {2., 0.}})};
        auto box = parameter_box(*cone, loops);
        ASSERT_NEAR(box.v1, apex, 1e-12);
        auto n = cone->to_nurbs(box); // degenerate side at the apex
        ASSERT_LT(dist(n->value(cone->nurbs_parameters(box, 1.5, apex)[0], apex), cone->value(0., apex)), 1e-12);
    }
    // a plane: extent plus 10 %
    {
        auto pl = r.surface(20);
        std::vector<std::vector<Point<3>>> loops{loop(*pl, {{0., 0.}, {2., 0.}, {2., 1.}, {0., 1.}})};
        auto box = parameter_box(*pl, loops);
        ASSERT_NEAR(box.u1, -0.2, 1e-12);
        ASSERT_NEAR(box.u2, 2.2, 1e-12);
        ASSERT_NEAR(box.v1, -0.1, 1e-12);
    }
}

TEST(tests_step_geometry, invalid_content)
{
    auto f = parse(step_file(frame + R"(
#20=CYLINDRICAL_SURFACE('',#13,-2.);
#21=CONICAL_SURFACE('',#13,1.,1.6);
#22=CIRCLE('',#13,0.);
#23=DIRECTION('',(0.,0.,0.));
#24=AXIS2_PLACEMENT_3D('',#10,#23,#12);
#25=PLANE('',#24);
#26=AXIS2_PLACEMENT_3D('',#10,#11,#11);
#27=PLANE('',#26);
#28=B_SPLINE_CURVE_WITH_KNOTS('',2,(#10,#10),.UNSPECIFIED.,.F.,.F.,(3,3),(0.,1.),.UNSPECIFIED.);
#29=FOO_SURFACE('',#13);
#30=PLANE('',#10);
)"));
    GeometryReader r{f, Units{}};
    ASSERT_THROW(r.surface(20), StepError);       // negative radius
    ASSERT_THROW(r.surface(21), StepError);       // semi angle beyond pi / 2
    ASSERT_THROW(r.curve<3>(22), StepError);      // null radius
    ASSERT_THROW(r.surface(25), StepError);       // null direction
    ASSERT_THROW(r.surface(27), StepError);       // reference along the axis
    ASSERT_THROW(r.curve<3>(28), StepError);      // knots do not match poles
    ASSERT_THROW(r.surface(29), StepUnsupported); // unknown surface
    ASSERT_THROW(r.surface(30), StepUnsupported); // a point where a placement is expected
    ASSERT_THROW(r.curve<3>(999), P21AccessError); // dangling reference
    ASSERT_THROW(r.surface(13), StepUnsupported);  // a placement where a surface is expected
}
