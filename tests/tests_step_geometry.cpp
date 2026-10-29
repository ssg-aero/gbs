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

    // Millimetres and degrees, the degree factor written with 10 digits as SolidWorks does:
    // 360 x 0.01745329252 = 2 pi + 2e-11.
    const std::string ctx_degrees_10_digits = R"(
#1=(GEOMETRIC_REPRESENTATION_CONTEXT(3)GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT((#5))GLOBAL_UNIT_ASSIGNED_CONTEXT((#2,#6,#4))REPRESENTATION_CONTEXT('',''));
#2=(LENGTH_UNIT()NAMED_UNIT(*)SI_UNIT(.MILLI.,.METRE.));
#3=(NAMED_UNIT(*)PLANE_ANGLE_UNIT()SI_UNIT($,.RADIAN.));
#4=(NAMED_UNIT(*)SI_UNIT($,.STERADIAN.)SOLID_ANGLE_UNIT());
#5=UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(1.E-05),#2,'distance_accuracy_value','');
#6=(CONVERSION_BASED_UNIT('DEGREE',#7)NAMED_UNIT(#8)PLANE_ANGLE_UNIT());
#7=PLANE_ANGLE_MEASURE_WITH_UNIT(PLANE_ANGLE_MEASURE(0.01745329252),#3);
#8=DIMENSIONAL_EXPONENTS(0.,0.,0.,0.,0.,0.,0.);
)";

    // True if f throws a StepError reporting a cyclic reference.
    template <typename F>
    bool throws_cycle(F &&f)
    {
        try
        {
            f();
        }
        catch (const StepError &e)
        {
            return std::string(e.what()).find("cyclic reference") != std::string::npos;
        }
        return false;
    }

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

TEST(tests_step_geometry, units_as_simple_instances)
{
    // simple SI_UNIT and CONVERSION_BASED_UNIT instances carry no LENGTH_UNIT nor PLANE_ANGLE_UNIT:
    // an SI unit is classified by its name, a conversion based unit by the unit of its factor
    {
        auto f = parse(step_file("", R"(
#1=(GEOMETRIC_REPRESENTATION_CONTEXT(3)GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT((#9))GLOBAL_UNIT_ASSIGNED_CONTEXT((#2,#3,#4,#5))REPRESENTATION_CONTEXT('',''));
#2=SI_UNIT(*,$,.METRE.);
#3=CONVERSION_BASED_UNIT(#8,'DEGREE',#7);
#4=SI_UNIT(*,$,.STERADIAN.);
#5=SI_UNIT(*,$,.SECOND.);
#6=SI_UNIT(*,$,.RADIAN.);
#7=PLANE_ANGLE_MEASURE_WITH_UNIT(PLANE_ANGLE_MEASURE(0.0174532925199433),#6);
#8=DIMENSIONAL_EXPONENTS(0.,0.,0.,0.,0.,0.,0.);
#9=UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(1.E-06),#2,'distance_accuracy_value','');
)"));
        auto u = read_units(f, 1);
        ASSERT_EQ(u.length_mm, 1000.);
        ASSERT_EQ(u.length, 1000.);
        ASSERT_EQ(u.length_name, "METRE");
        ASSERT_NEAR(u.angle, pi / 180., 1e-15);
        ASSERT_EQ(u.angle_name, "DEGREE");
        ASSERT_NEAR(*u.uncertainty, 1e-3, 1e-15); // 1e-6 m
    }
    {
        auto f = parse(step_file("", R"(
#1=(GEOMETRIC_REPRESENTATION_CONTEXT(3)GLOBAL_UNIT_ASSIGNED_CONTEXT((#2,#3))REPRESENTATION_CONTEXT('',''));
#2=SI_UNIT(*,.CENTI.,.METRE.);
#3=SI_UNIT(*,$,.RADIAN.);
)"));
        auto u = read_units(f, 1);
        ASSERT_NEAR(u.length_mm, 10., 1e-12);
        ASSERT_EQ(u.length_name, "CENTI METRE");
        ASSERT_EQ(u.angle, 1.);
        ASSERT_EQ(u.angle_name, "RADIAN");
    }
}

TEST(tests_step_geometry, axis2_placement_defaults)
{
    // ISO 10303-42 build_axes / first_proj_axis: axis (0, 0, 1) by default; reference (1, 0, 0) by
    // default, (0, 1, 0) when the axis is (1, 0, 0) or (-1, 0, 0); y = z x x
    auto f = parse(step_file(R"(
#10=CARTESIAN_POINT('',(0.,0.,0.));
#11=DIRECTION('',(1.,0.,0.));
#12=AXIS2_PLACEMENT_3D('',#10,#11,$);
#13=DIRECTION('',(-2.,0.,0.));
#14=AXIS2_PLACEMENT_3D('',#10,#13,$);
#15=DIRECTION('',(0.,1.,0.));
#16=AXIS2_PLACEMENT_3D('',#10,#15,$);
#17=AXIS2_PLACEMENT_3D('',#10,$,$);
#18=CARTESIAN_POINT('',(0.,0.));
#19=AXIS2_PLACEMENT_2D('',#18,$);
#20=AXIS1_PLACEMENT('',#10,$);
#21=CIRCLE('',#12,10.);
)"));
    GeometryReader r{f, Units{}};
    const Point<3> X{1., 0., 0.}, Y{0., 1., 0.}, Z{0., 0., 1.};
    const auto f1 = r.frame(12);
    ASSERT_LT(dist(f1.x, Y), tol);
    ASSERT_LT(dist(f1.y, Z), tol);
    const auto f2 = r.frame(14);
    ASSERT_LT(dist(f2.z, -1. * X), tol);
    ASSERT_LT(dist(f2.x, Y), tol);
    ASSERT_LT(dist(f2.y, -1. * Z), tol);
    const auto f3 = r.frame(16);
    ASSERT_LT(dist(f3.x, X), tol);
    ASSERT_LT(dist(f3.y, -1. * Z), tol);
    const auto f4 = r.frame(17);
    ASSERT_LT(dist(f4.x, X), tol);
    ASSERT_LT(dist(f4.y, Y), tol);
    ASSERT_LT(dist(f4.z, Z), tol);
    const auto [o2, x2, y2] = r.frame2d(19);
    ASSERT_EQ(x2, (Point<2>{1., 0.}));
    ASSERT_EQ(y2, (Point<2>{0., 1.}));
    ASSERT_LT(dist(r.axis1(20).second, Z), tol);
    // a circle on the frame of axis (1, 0, 0): starts on y, turns towards z
    auto c = r.curve<3>(21);
    ASSERT_LT(dist(c->value(0.), Point<3>{0., 10., 0.}), 1e-12);
    ASSERT_LT(dist(c->value(pi / 2), Point<3>{0., 0., 10.}), 1e-12);
}

TEST(tests_step_geometry, full_turn_with_a_rounded_degree_factor)
{
    // 360 degrees are 2 pi + 2e-11 rad: a full turn, neither an arc of 2e-11 rad nor a span above 2 pi
    auto f = parse(step_file(frame + R"(
#20=CIRCLE('',#13,10.);
#21=TRIMMED_CURVE('',#20,(PARAMETER_VALUE(0.)),(PARAMETER_VALUE(360.)),.T.,.PARAMETER.);
#22=TRIMMED_CURVE('',#20,(PARAMETER_VALUE(90.)),(PARAMETER_VALUE(450.)),.T.,.PARAMETER.);
#23=TRIMMED_CURVE('',#20,(PARAMETER_VALUE(360.)),(PARAMETER_VALUE(0.)),.F.,.PARAMETER.);
#24=CYLINDRICAL_SURFACE('',#13,10.);
#25=RECTANGULAR_TRIMMED_SURFACE('',#24,0.,360.,0.,5.,.T.,.T.);
#26=TOROIDAL_SURFACE('',#13,20.,5.);
#27=RECTANGULAR_TRIMMED_SURFACE('',#26,90.,450.,-180.,180.,.T.,.T.);
#28=TRIMMED_CURVE('',#20,(PARAMETER_VALUE(0.)),(PARAMETER_VALUE(359.999)),.T.,.PARAMETER.);
#29=TRIMMED_CURVE('',#20,(PARAMETER_VALUE(0.)),(PARAMETER_VALUE(0.001)),.T.,.PARAMETER.);
)",
                             ctx_degrees_10_digits));
    GeometryReader r{f, read_units(f, 1)};
    ASSERT_GT(360. * r.units().angle, 2 * pi + 1e-11);
    const double turn = 2 * pi * (1. + 1e-15); // a turn up to rounding
    for (std::uint64_t id : {21, 22, 23})
    {
        auto c = r.curve<3>(id);
        const auto [t1, t2] = c->domain();
        ASSERT_NEAR(t2 - t1, 2 * pi, 1e-14);
        ASSERT_LE(t2 - t1, turn);
        auto n = c->to_nurbs(t1, t2);
        ASSERT_LT(dist(n.curve->begin(), n.curve->end()), 1e-12);
        ASSERT_LT(dist(n.curve->value(c->nurbs_parameter(t1, t2, t1 + pi)), c->value(t1 + pi)), 1e-12);
    }
    ASSERT_TRUE(r.curve<3>(23)->reversed());
    // trimmed cylinder and torus: a turn in u (and in v for the torus); NURBS and primitive build
    for (std::uint64_t id : {25, 27})
    {
        auto s = r.surface(id);
        const auto d = s->domain();
        ASSERT_NEAR(d.u2 - d.u1, 2 * pi, 1e-14);
        ASSERT_LE(d.u2 - d.u1, turn);
        auto n = s->to_nurbs(d);
        auto e = s->to_surface(d);
        const double u = d.u1 + 1., v = d.v1 + 0.5;
        auto [nu, nv] = s->nurbs_parameters(d, u, v);
        ASSERT_LT(dist(n->value(nu, nv), s->value(u, v)), 1e-11);
        ASSERT_LT(dist(e->value(u, v), s->value(u, v)), 1e-11);
    }
    const auto dt = r.surface(27)->domain();
    ASSERT_NEAR(dt.v2 - dt.v1, 2 * pi, 1e-14);
    ASSERT_LE(dt.v2 - dt.v1, turn);
    // an arc 0.001 degree short of a turn, or 0.001 degree long, stays what it is
    const auto a = r.curve<3>(28)->domain();
    ASSERT_NEAR(a[1] - a[0], 359.999 * r.units().angle, 1e-12);
    const auto b = r.curve<3>(29)->domain();
    ASSERT_NEAR(b[1] - b[0], 0.001 * r.units().angle, 1e-15);
}

TEST(tests_step_geometry, rectangular_trimmed_surface_senses)
{
    // ISO 10303-42: in a periodic direction the range runs from the first bound to the second in the
    // given sense, across the seam if need be; in another direction the bounds are only sorted
    auto f = parse(step_file(frame + R"(
#20=CYLINDRICAL_SURFACE('',#13,10.);
#21=RECTANGULAR_TRIMMED_SURFACE('',#20,350.,10.,0.,5.,.T.,.T.);
#22=RECTANGULAR_TRIMMED_SURFACE('',#20,10.,350.,5.,0.,.F.,.F.);
#23=RECTANGULAR_TRIMMED_SURFACE('',#20,10.,350.,0.,5.,.T.,.T.);
#24=RECTANGULAR_TRIMMED_SURFACE('',#20,350.,10.,0.,5.,.F.,.T.);
#30=TOROIDAL_SURFACE('',#13,20.,5.);
#31=RECTANGULAR_TRIMMED_SURFACE('',#30,0.,90.,300.,60.,.T.,.T.);
#32=RECTANGULAR_TRIMMED_SURFACE('',#30,0.,90.,60.,300.,.T.,.F.);
#40=PLANE('',#13);
#41=RECTANGULAR_TRIMMED_SURFACE('',#40,3.,-1.,2.,1.,.F.,.F.);
)"));
    GeometryReader r{f, Units{.angle = pi / 180.}};
    const double deg = pi / 180.;
    const std::vector<std::pair<std::uint64_t, UVBox>> cases{
        {21, {350. * deg, 370. * deg, 0., 5.}},  // forward across the seam: 20 degrees
        {22, {350. * deg, 370. * deg, 0., 5.}},  // backward from 10 to 350: the same 20 degrees
        {23, {10. * deg, 350. * deg, 0., 5.}},   // forward from 10 to 350: 340 degrees
        {24, {10. * deg, 350. * deg, 0., 5.}},   // backward from 350 to 10
        {31, {0., 90. * deg, 300. * deg, 420. * deg}},
        {32, {0., 90. * deg, 300. * deg, 420. * deg}},
        {41, {-1., 3., 1., 2.}}, // plane: sorted
    };
    for (const auto &[id, b] : cases)
    {
        auto s = r.surface(id);
        const auto d = s->domain();
        ASSERT_NEAR(d.u1, b.u1, 1e-12);
        ASSERT_NEAR(d.u2, b.u2, 1e-12);
        ASSERT_NEAR(d.v1, b.v1, 1e-12);
        ASSERT_NEAR(d.v2, b.v2, 1e-12);
        // exact NURBS and primitive on the range; points inverted into the range
        auto n = s->to_nurbs(d);
        auto e = s->to_surface(d);
        for (double i : {0.1, 0.5, 0.9})
        {
            const double u = d.u1 + i * (d.u2 - d.u1), v = d.v1 + (1. - i) * (d.v2 - d.v1);
            const auto p = s->value(u, v);
            auto [nu, nv] = s->nurbs_parameters(d, u, v);
            ASSERT_LT(dist(n->value(nu, nv), p), 1e-11);
            ASSERT_LT(dist(e->value(u, v), p), 1e-11);
            const auto uv = *s->parameters(p);
            ASSERT_NEAR(uv[0], u, 1e-11);
            ASSERT_NEAR(uv[1], v, 1e-11);
        }
    }
    // the 20 degree patch is centred on the reference direction
    ASSERT_LT(dist(r.surface(21)->value(0., 2.), Point<3>{11., 2., 5.}), 1e-12);
}

TEST(tests_step_geometry, bspline_inversion_hairpin_and_ties)
{
    auto f = parse(step_file(R"(
#10=CARTESIAN_POINT('',(0.,0.,0.));
#11=CARTESIAN_POINT('',(100.,0.,0.));
#12=CARTESIAN_POINT('',(100.,0.005,0.));
#13=CARTESIAN_POINT('',(0.,0.005,0.));
#20=B_SPLINE_CURVE_WITH_KNOTS('',1,(#10,#11,#12,#13),.UNSPECIFIED.,.F.,.F.,(2,1,1,2),(0.,0.3,0.32,1.),.UNSPECIFIED.);
#30=CARTESIAN_POINT('',(10.,0.,0.));
#31=CARTESIAN_POINT('',(10.,10.,0.));
#32=CARTESIAN_POINT('',(-10.,10.,0.));
#33=CARTESIAN_POINT('',(-10.,-10.,0.));
#34=CARTESIAN_POINT('',(10.,-10.,0.));
#35=CARTESIAN_POINT('',(10.,0.,0.));
#40=B_SPLINE_CURVE_WITH_KNOTS('',3,(#30,#31,#32,#33,#34,#35),.UNSPECIFIED.,.T.,.F.,(4,2,4),(0.,0.5,1.),.UNSPECIFIED.);
#50=CARTESIAN_POINT('',(0.,0.,0.));
#51=CARTESIAN_POINT('',(1.,3.,0.));
#52=CARTESIAN_POINT('',(2.,-3.,1.));
#53=CARTESIAN_POINT('',(3.,3.,0.));
#54=CARTESIAN_POINT('',(4.,-3.,1.));
#55=CARTESIAN_POINT('',(5.,0.,0.));
#56=BEZIER_CURVE('',5,(#50,#51,#52,#53,#54,#55),.UNSPECIFIED.,.F.,.F.);
)"));
    GeometryReader r{f, Units{}};
    // a degree 1 curve going to x = 100 and coming back 0.005 mm higher: each branch has its samples
    auto c = r.curve<3>(20);
    for (int i = 0; i <= 400; ++i)
    {
        const double t = i / 400.;
        ASSERT_NEAR(c->parameter(c->value(t)), t, 1e-12);
    }
    // halfway between the branches, equally near both: the smaller parameter
    ASSERT_LT(c->parameter(Point<3>{50., 0.0025, 0.}), 0.3);
    // a closed curve: its closure point gives the start of the domain
    auto closed = r.curve<3>(40);
    ASSERT_NEAR(closed->parameter(Point<3>{10., 0., 0.}), 0., 1e-12);
    for (int i = 1; i < 100; ++i)
        ASSERT_NEAR(closed->parameter(closed->value(i / 100.)), i / 100., 1e-12);
    // a wavy quintic Bezier (one knot span), and points off the curve: the nearest point
    auto w = r.curve<3>(56);
    for (int i = 0; i <= 200; ++i)
    {
        const double t = i / 200.;
        ASSERT_NEAR(w->parameter(w->value(t)), t, 1e-11);
        const auto off = w->value(t) + 0.01 * Point<3>{0., 0., 1.};
        const double s = w->parameter(off);
        ASSERT_LE(dist(w->value(s), off), dist(w->value(t), off) + 1e-12);
    }
}

TEST(tests_step_geometry, extrusion_inversion_with_an_oblique_vector)
{
    // S(u, v) = C(u) + v V with V not normal to the plane of C: u minimises the distance seen along V
    auto f = parse(step_file(frame + R"(
#20=CARTESIAN_POINT('',(0.,0.,0.));
#21=CARTESIAN_POINT('',(5.,5.,5.));
#22=CARTESIAN_POINT('',(10.,0.,10.));
#23=B_SPLINE_CURVE_WITH_KNOTS('',2,(#20,#21,#22),.UNSPECIFIED.,.F.,.F.,(3,3),(0.,1.),.UNSPECIFIED.);
#24=DIRECTION('',(0.,0.,1.));
#25=VECTOR('',#24,1.);
#26=SURFACE_OF_LINEAR_EXTRUSION('',#23,#25);
#30=CIRCLE('',#13,2.);
#31=TRIMMED_CURVE('',#30,(PARAMETER_VALUE(300.)),(PARAMETER_VALUE(60.)),.T.,.PARAMETER.);
#32=DIRECTION('',(1.,1.,1.));
#33=VECTOR('',#32,2.);
#34=SURFACE_OF_LINEAR_EXTRUSION('',#31,#33);
#35=TRIMMED_CURVE('',#23,(PARAMETER_VALUE(0.2)),(PARAMETER_VALUE(0.9)),.T.,.PARAMETER.);
#36=SURFACE_OF_LINEAR_EXTRUSION('',#35,#33);
)"));
    GeometryReader r{f, Units{.angle = pi / 180.}};
    for (std::uint64_t id : {26, 34, 36})
    {
        auto s = r.surface(id);
        const auto d = s->domain();
        for (int i = 0; i <= 10; ++i)
            for (double v : {-3., 0., 2., 4.})
            {
                const double u = d.u1 + (d.u2 - d.u1) * i / 10.;
                const auto p = s->value(u, v);
                const auto uv = *s->parameters(p);
                ASSERT_NEAR(uv[0], u, 1e-9);
                ASSERT_NEAR(uv[1], v, 1e-9);
            }
    }
}

TEST(tests_step_geometry, revolution_inversion_skew_generatrix)
{
    // v minimises the distance in the meridian half-plane, u turns the generatrix point onto p
    auto f = parse(step_file(R"(
#10=CARTESIAN_POINT('',(0.,0.,0.));
#11=DIRECTION('',(0.,0.,1.));
#12=AXIS1_PLACEMENT('',#10,#11);
#20=CARTESIAN_POINT('',(5.,0.,0.));
#21=DIRECTION('',(0.,1.,1.));
#22=VECTOR('',#21,1.41421356237310);
#23=LINE('',#20,#22);
#24=SURFACE_OF_REVOLUTION('',#23,#12);
#30=CARTESIAN_POINT('',(-5.,0.,0.));
#31=DIRECTION('',(1.,0.,1.));
#32=VECTOR('',#31,1.41421356237310);
#33=LINE('',#30,#32);
#34=SURFACE_OF_REVOLUTION('',#33,#12);
#40=CARTESIAN_POINT('',(4.,0.,0.));
#41=CARTESIAN_POINT('',(6.,3.,2.));
#42=CARTESIAN_POINT('',(3.,-2.,5.));
#43=CARTESIAN_POINT('',(5.,1.,8.));
#44=B_SPLINE_CURVE_WITH_KNOTS('',3,(#40,#41,#42,#43),.UNSPECIFIED.,.F.,.F.,(4,4),(0.,1.),.UNSPECIFIED.);
#45=SURFACE_OF_REVOLUTION('',#44,#12);
#50=CARTESIAN_POINT('',(1.,0.,0.));
#51=DIRECTION('',(0.,1.,0.));
#52=DIRECTION('',(1.,0.,0.));
#53=AXIS2_PLACEMENT_3D('',#50,#51,#52);
#54=CIRCLE('',#53,3.);
#55=SURFACE_OF_REVOLUTION('',#54,#12);
#56=TRIMMED_CURVE('',#33,(PARAMETER_VALUE(2.)),(PARAMETER_VALUE(9.)),.T.,.PARAMETER.);
#57=SURFACE_OF_REVOLUTION('',#56,#12);
)"));
    GeometryReader r{f, Units{}};
    // hyperboloid of one sheet (skew line), double cone (line through the axis at v = 5), skew cubic,
    // circle across the axis, trimmed line through the axis
    for (auto [id, v1, v2] : std::vector<std::tuple<std::uint64_t, double, double>>{
             {24, -3., 6.}, {34, -2., 12.}, {45, 0., 1.}, {55, 0., 2 * pi}, {57, 2., 9.}})
    {
        auto s = r.surface(id);
        for (int i = 0; i < 8; ++i)
            for (int j = 0; j <= 12; ++j)
            {
                const double u = 2 * pi * (i + 0.3) / 8., v = v1 + (v2 - v1) * j / 12.;
                const auto p = s->value(u, v);
                const auto uv = *s->parameters(p);
                ASSERT_LT(dist(s->value(uv[0], uv[1]), p), 1e-9);
                if ((id == 24 || id == 34 || id == 45) && std::hypot(p[0], p[1]) > 1e-9) // unique parameters
                {
                    ASSERT_NEAR(uv[0], u, 1e-9);
                    ASSERT_NEAR(uv[1], v, 1e-9);
                }
            }
    }
    // a face of the hyperboloid, u in [0, 1], v in [0, 3]: its box contains it
    auto s = r.surface(24);
    std::vector<Point<3>> loop;
    for (int i = 0; i < 16; ++i)
        loop.push_back(s->value(i / 16., 0.));
    for (int i = 0; i < 16; ++i)
        loop.push_back(s->value(1., 3. * i / 16.));
    for (int i = 0; i < 16; ++i)
        loop.push_back(s->value(1. - i / 16., 3.));
    for (int i = 0; i < 16; ++i)
        loop.push_back(s->value(0., 3. - 3. * i / 16.));
    std::vector<std::vector<Point<3>>> loops{loop};
    const auto b = parameter_box(*s, loops);
    ASSERT_LE(b.u1, 0.);
    ASSERT_GE(b.u2, 1.);
    ASSERT_LE(b.v1, 0.);
    ASSERT_GE(b.v2, 3.);
}

TEST(tests_step_geometry, offset_cone_keeps_the_parametrization)
{
    // ISO: S'(u, v) = S(u, v) + d N(u, v), N = cos a X(u) - sin a Z for the cone
    auto f = parse(step_file(frame + R"(
#20=CONICAL_SURFACE('',#13,10.,30.);
#21=OFFSET_SURFACE('',#20,2.,.F.);
#22=OFFSET_SURFACE('',#20,-3.,.F.);
#23=CONICAL_SURFACE('',#13,10.,-20.);
#24=OFFSET_SURFACE('',#23,1.5,.F.);
#25=RECTANGULAR_TRIMMED_SURFACE('',#21,0.,90.,0.,10.,.T.,.T.);
#26=OFFSET_SURFACE('',#20,-20.,.F.);
)"));
    const double deg = pi / 180.;
    GeometryReader r{f, Units{.angle = deg}};
    for (auto [bid, oid, d, a] : std::vector<std::tuple<std::uint64_t, std::uint64_t, double, double>>{
             {20, 21, 2., 30. * deg}, {20, 22, -3., 30. * deg}, {23, 24, 1.5, -20. * deg}})
    {
        auto base = r.surface(bid), off = r.surface(oid);
        for (double u : {0., 0.4, 2., 5.})
            for (double v : {-5., 0., 5., 12.})
            {
                const Point<3> n = std::cos(a) * Point<3>{std::cos(u), std::sin(u), 0.} - std::sin(a) * Point<3>{0., 0., 1.};
                ASSERT_LT(dist(off->value(u, v), base->value(u, v) + d * n), 1e-12);
            }
    }
    // trimmed by parameters on the offset: the edge v = 0 is at z = 3 - d sin a
    auto t = r.surface(25);
    const auto dm = t->domain();
    ASSERT_NEAR(t->value(dm.u1, dm.v1)[2], 3. - 2. * std::sin(30. * deg), 1e-12);
    // inwards beyond the reference circle (R + d cos a < 0): not a cone of the same parameters
    ASSERT_THROW(r.surface(26), StepUnsupported);
}

TEST(tests_step_geometry, degenerate_torus_apple_and_lemon)
{
    // minor radius 5 > major radius 2: select_outer .T. is the apple (v in [-t, t]), .F. the lemon
    // (v in [t, 2 pi - t], where R + r cos v < 0), cos t = -R / r
    auto f = parse(step_file(frame + R"(
#20=DEGENERATE_TOROIDAL_SURFACE('',#13,2.,5.,.T.);
#21=DEGENERATE_TOROIDAL_SURFACE('',#13,2.,5.,.F.);
#22=OFFSET_SURFACE('',#21,0.5,.F.);
)"));
    GeometryReader r{f, Units{}};
    const double t = std::acos(-2. / 5.);
    const std::vector<std::pair<std::uint64_t, std::array<double, 2>>> parts{{20, {-t, t}}, {21, {t, 2 * pi - t}}};
    for (const auto &[id, range] : parts)
    {
        auto s = r.surface(id);
        const auto d = s->domain();
        ASSERT_EQ(d.u1, 0.);
        ASSERT_EQ(d.u2, 2 * pi);
        ASSERT_NEAR(d.v1, range[0], 1e-15);
        ASSERT_NEAR(d.v2, range[1], 1e-15);
        ASSERT_TRUE(s->periodic()[0]);
        ASSERT_FALSE(s->periodic()[1]);
        auto n = s->to_nurbs(d);
        auto e = s->to_surface(d);
        for (int i = 0; i < 8; ++i)
            for (int j = 0; j <= 10; ++j)
            {
                const double u = 2 * pi * (i + 0.3) / 8., v = d.v1 + (d.v2 - d.v1) * j / 10.;
                const auto p = s->value(u, v);
                auto [nu, nv] = s->nurbs_parameters(d, u, v);
                ASSERT_LT(dist(n->value(nu, nv), p), 1e-11);
                ASSERT_LT(dist(e->value(u, v), p), 1e-11);
                const auto uv = *s->parameters(p);
                ASSERT_LT(dist(s->value(uv[0], uv[1]), p), 1e-11);
                ASSERT_GE(uv[1], d.v1 - 1e-12);
                ASSERT_LE(uv[1], d.v2 + 1e-12);
                if (j != 0 && j != 10) // the ends are the poles, on the axis
                {
                    ASSERT_NEAR(uv[0], u, 1e-11);
                    ASSERT_NEAR(uv[1], v, 1e-11);
                }
            }
    }
    // the offset of a part changes its range of v: not supported rather than wrong
    ASSERT_THROW(r.surface(22), StepUnsupported);
}

TEST(tests_step_geometry, cyclic_references)
{
    auto f = parse(step_file(frame + R"(
#20=TRIMMED_CURVE('',#20,(PARAMETER_VALUE(0.)),(PARAMETER_VALUE(1.)),.T.,.PARAMETER.);
#21=SURFACE_CURVE('',#21,(),.CURVE_3D.);
#22=COMPOSITE_CURVE_SEGMENT(.CONTINUOUS.,.T.,#23);
#23=COMPOSITE_CURVE('',(#22),.F.);
#24=TRIMMED_CURVE('',#25,(PARAMETER_VALUE(0.)),(PARAMETER_VALUE(1.)),.T.,.PARAMETER.);
#25=SURFACE_CURVE('',#24,(),.CURVE_3D.);
#30=RECTANGULAR_TRIMMED_SURFACE('',#30,0.,1.,0.,1.,.T.,.T.);
#31=OFFSET_SURFACE('',#32,1.,.F.);
#32=OFFSET_SURFACE('',#31,1.,.F.);
#40=CIRCLE('',#13,1.);
)"));
    GeometryReader r{f, Units{}};
    for (std::uint64_t id : {20, 21, 23, 24, 25})
    {
        ASSERT_TRUE(throws_cycle([&] { (void)r.curve<3>(id); }));
        ASSERT_TRUE(throws_cycle([&] { (void)r.curve<2>(id); }));
    }
    for (std::uint64_t id : {30, 31, 32})
        ASSERT_TRUE(throws_cycle([&] { (void)r.surface(id); }));
    // the reader stays usable, and reports the cycle again
    ASSERT_LT(dist(r.curve<3>(40)->value(0.), Point<3>{2., 2., 3.}), tol);
    ASSERT_TRUE(throws_cycle([&] { (void)r.curve<3>(20); }));
}

TEST(tests_step_geometry, bsplines_with_too_few_poles)
{
    // a degree of 1 at least and degree + 1 poles at least are checked before any knot is generated
    auto f = parse(step_file(R"(
#10=CARTESIAN_POINT('',(0.,0.,0.));
#11=CARTESIAN_POINT('',(1.,0.,0.));
#12=CARTESIAN_POINT('',(1.,1.,0.));
#20=QUASI_UNIFORM_CURVE('',3,(#10,#11),.UNSPECIFIED.,.F.,.F.);
#21=UNIFORM_CURVE('',3,(#10,#11),.UNSPECIFIED.,.F.,.F.);
#22=BEZIER_CURVE('',2,(),.UNSPECIFIED.,.F.,.F.);
#23=B_SPLINE_CURVE_WITH_KNOTS('',3,(#10,#11),.UNSPECIFIED.,.F.,.F.,(3,3),(0.,1.),.UNSPECIFIED.);
#24=QUASI_UNIFORM_CURVE('',0,(#10,#11),.UNSPECIFIED.,.F.,.F.);
#25=QUASI_UNIFORM_CURVE('',-1,(#10,#11),.UNSPECIFIED.,.F.,.F.);
#26=B_SPLINE_CURVE_WITH_KNOTS('',1,(#10,#11,#12),.UNSPECIFIED.,.F.,.F.,(2,-1,2),(0.,0.5,1.),.UNSPECIFIED.);
#27=B_SPLINE_CURVE_WITH_KNOTS('',1,(#10,#11,#12),.UNSPECIFIED.,.F.,.F.,(2,1000000000,2),(0.,0.5,1.),.UNSPECIFIED.);
#30=QUASI_UNIFORM_SURFACE('',1,3,((#10,#11),(#11,#12)),.UNSPECIFIED.,.F.,.F.,.F.);
#31=BEZIER_SURFACE('',3,1,((#10,#11),(#11,#12)),.UNSPECIFIED.,.F.,.F.,.F.);
#32=B_SPLINE_SURFACE_WITH_KNOTS('',1,3,((#10,#11),(#11,#12)),.UNSPECIFIED.,.F.,.F.,.F.,(2,2),(3,3),(0.,1.),(0.,1.),.UNSPECIFIED.);
#40=QUASI_UNIFORM_CURVE('',1,(#10,#11),.UNSPECIFIED.,.F.,.F.);
)"));
    GeometryReader r{f, Units{}};
    for (std::uint64_t id : {20, 21, 22, 23, 24, 25, 26, 27})
        ASSERT_THROW(r.curve<3>(id), StepError);
    for (std::uint64_t id : {30, 31, 32})
        ASSERT_THROW(r.surface(id), StepError);
    // the smallest valid one: a segment
    ASSERT_LT(dist(r.curve<3>(40)->value(0.5), Point<3>{0.5, 0., 0.}), tol);
}
