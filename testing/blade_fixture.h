#pragma once
// A blade on its hub, for the tests of the fillets (gbs and against OCCT): the foil of
// tests_foils (type2_blade: a camber line thickened by a NACA law, its leading and trailing
// edges rounded by C3 connections), stacked and twisted along z from below the hub to its tip,
// skinned and closed by two planes; the hub a cylinder along x of radius 1. The blade fused to
// the hub: its root a closed edge round the blade on the cylinder, the edge to fillet.
#include <gbs-brep/brep>
#include <gbs/bscapprox.h>
#include <gbs/bscbuild.h>
#include <gbs/bscanalysis.h>
#include <gbs/bsctools.h>
#include <gbs/bssbuild.h>

#include <cmath>
#include <memory>
#include <numbers>
#include <vector>

namespace blade_fixture
{
    using T = double;
    using namespace gbs;
    using namespace gbs::brep;

    /// The closed foil, chord about 1 along x from (0, 0): one B-spline, its start and end at the trailing edge.
    inline auto foil() -> BSCurve<T, 2>
    {
        const T b1 = std::numbers::pi / 24., t1 = 0.3, t2 = 0.3;
        const BSCurve<T, 2> camber{points_vector<T, 2>{{0., 0.}, {t1 * std::cos(b1), t1 * std::sin(b1)}, {1. - t2 * std::cos(b1), t2 * std::sin(b1)}, {1., 0.}},
                                   std::vector<T>{0., 0., 0., 0.5, 1., 1., 1.}, 2};
        const T chord = length(camber, 0., 1.);
        auto naca = [](T x) { return 5. * 0.1 * (0.2969 * std::sqrt(x) - 0.1260 * x - 0.3516 * x * x + 0.2843 * x * x * x - 0.1015 * x * x * x * x); };
        std::vector<T> u(100);
        for (std::size_t i = 0; i < u.size(); ++i)
            u[i] = T(i) / T(u.size() - 1);
        auto side = [&](T sign, std::size_t n_poles) {
            points_vector<T, 2> pts;
            for (T x : u)
            {
                auto t = camber.value(x, 1);
                const point<T, 2> n{-t[1] / norm(t), t[0] / norm(t)};
                const T l = x > 0. ? length(camber, 0., x) / chord : 0.;
                const T ep = naca(0.01 * (1. - l) + 0.95 * l);
                pts.push_back(camber.value(x) + (sign * ep) * n);
            }
            return approx(pts, 5, n_poles, u, true);
        };
        auto s1 = side(1., 10), s2 = side(-1., 12);
        s1.reverse();
        const auto le = c3_connect(s1, s2), te = c3_connect(s2, s1);
        return join(join(join(s1, le), s2), te);
    }

    /// The section of the foil at height z, turned by `twist` about (0.5, 0).
    inline auto section(const BSCurve<T, 2> &f, T z, T twist) -> std::shared_ptr<BSCurve<T, 3>>
    {
        const T c = std::cos(twist), s = std::sin(twist);
        points_vector<T, 3> poles;
        for (const auto &p : f.poles())
            poles.push_back({0.5 + c * (p[0] - 0.5) - s * p[1], s * (p[0] - 0.5) + c * p[1], z});
        return std::make_shared<BSCurve<T, 3>>(poles, f.knotsFlats(), f.degree());
    }

    /// The solid blade through the sections at the heights zs, turned by `twists`, closed by planes.
    inline auto blade(Model<T> &m, const std::vector<T> &zs, const std::vector<T> &twists) -> SolidId
    {
        const auto f = foil();
        std::vector<std::shared_ptr<Curve<T, 3>>> secs;
        for (std::size_t i = 0; i < zs.size(); ++i)
            secs.push_back(section(f, zs[i], twists[i]));
        const auto skin = std::make_shared<BSSurface<T, 3>>(loft<T, 3>(secs, 3));
        std::vector<FaceId> faces{unwrap(make_face(m, std::shared_ptr<Surface<T, 3>>{skin}))};
        for (std::size_t i : {std::size_t{0}, zs.size() - 1})
        {
            const auto &c = secs[i];
            const auto [a, b] = c->bounds();
            const auto v = unwrap(make_vertex(m, c->value(a)));
            const auto e = unwrap(make_edge(m, c, a, b, v, v));
            const auto plane = std::make_shared<Plane<T>>(ax2<T, 3>{{{0., 0., zs[i]}, {0., 0., 1.}, {1., 0., 0.}}}, std::array<T, 4>{-5., 5., -5., 5.});
            faces.push_back(unwrap(make_face(m, plane, unwrap(make_wire(m, std::vector<EdgeId>{e})))));
        }
        const auto sewn = unwrap(sew(m, faces));
        return unwrap(make_solid(m, sewn.shells.front()));
    }

    /// The hub: the cylinder of radius 1 along x from x = -0.5 to 1.5.
    inline auto hub(Model<T> &m) -> SolidId
    {
        const ax2<T, 3> ax{{{-0.5, 0., 0.}, {1., 0., 0.}, {0., 1., 0.}}};
        const auto e = unwrap(make_edge(m, std::make_shared<EllipseCurve<T>>(1., 1., ax)));
        const auto f = unwrap(make_face(m, std::make_shared<Plane<T>>(ax, std::array<T, 4>{-2., 2., -2., 2.}), unwrap(make_wire(m, std::vector<EdgeId>{e}))));
        return explore<SolidId>(m, unwrap(make_prism(m, f, point<T, 3>{2., 0., 0.})).shape).front();
    }

    /// The blade fused to the hub, and the edges of its root (on the hub, round the blade).
    inline auto blade_on_hub(Model<T> &m) -> std::pair<SolidId, std::vector<EdgeId>>
    {
        const auto b = blade(m, {0.85, 1.2, 1.5, 1.8}, {-0.15, 0., 0.15, 0.3});
        const auto s = unwrap(fuse(m, ShapeId{hub(m)}, ShapeId{b})).solids.front();
        std::vector<EdgeId> root;
        for (auto e : explore<EdgeId>(m, ShapeId{s}))
        {
            const auto &ed = m.edge(e);
            if (ed.degenerate || !ed.curve)
                continue;
            const auto p = ed.curve->value((ed.u1 + ed.u2) / 2.);
            if (std::abs(std::hypot(p[1], p[2]) - 1.) < 1e-5 && std::abs(p[1]) < 0.5 && p[2] > 0. && p[0] > -0.3 && p[0] < 1.3)
                root.push_back(e);
        }
        return {s, root};
    }
}
