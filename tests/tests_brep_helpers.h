#pragma once
// Shared fixtures of the gbs-brep tests.

#include <gbs-brep/brep>
#include <gbs/bscbuild.h>

#include <array>
#include <memory>
#include <stdexcept>
#include <vector>

namespace brep_tests
{
    using T = double;
    using namespace gbs;
    using namespace gbs::brep;

    // A unit box built by hand, the way a STEP reader or a builder would:
    // 8 vertices, 12 edges (degree-1 B-spline segments), 6 planar degree-1
    // B-spline faces with their pcurves, 1 shell, 1 solid. Every face is
    // oriented with its normal pointing outwards, so the shell is closed,
    // manifold and consistently oriented with Forward face uses only.
    struct Box
    {
        Model<T> m;
        std::array<VertexId, 8> v;
        std::vector<EdgeId> e;
        std::array<FaceId, 6> f;
        ShellId shell;
        SolidId solid;
    };

    // Corner indices: bit 0 -> x, bit 1 -> y, bit 2 -> z
    inline point<T, 3> corner(unsigned i)
    {
        return {T(i & 1), T((i >> 1) & 1), T((i >> 2) & 1)};
    }

    inline EdgeId find_edge(const Model<T> &m, const std::vector<EdgeId> &edges, VertexId a, VertexId b, Orientation &orient)
    {
        for (auto eid : edges)
        {
            const auto &ed = m.edge(eid);
            if (ed.v1 == a && ed.v2 == b) { orient = Orientation::Forward; return eid; }
            if (ed.v1 == b && ed.v2 == a) { orient = Orientation::Reversed; return eid; }
        }
        throw std::runtime_error("edge not found");
    }

    inline Box make_box()
    {
        Box b;
        auto &m = b.m;

        for (unsigned i = 0; i < 8; ++i)
            b.v[i] = m.add(Vertex<T>{corner(i), 1e-7});

        // 12 edges along the axes, from the lower to the upper corner
        for (unsigned i = 0; i < 8; ++i)
            for (unsigned axis = 0; axis < 3; ++axis)
                if (!(i & (1u << axis)))
                {
                    unsigned j = i | (1u << axis);
                    auto crv = std::make_shared<BSCurve<T, 3>>(build_segment(corner(i), corner(j), true));
                    b.e.push_back(m.add(Edge<T>{crv, 0., 1., b.v[i], b.v[j], 1e-7}));
                }

        // 6 faces: for each axis, the two planes at coordinate 0 and 1.
        // Face corners in (u,v) order: (0,0) (1,0) (1,1) (0,1), surface S(u,v) = c00 + u (c10 - c00) + v (c01 - c00)
        size_t fi = 0;
        for (unsigned axis = 0; axis < 3; ++axis)
            for (unsigned side = 0; side < 2; ++side)
            {
                // (u, v) axes chosen so that the normal du ^ dv points outwards on both sides
                unsigned a1 = side ? (axis + 1) % 3 : (axis + 2) % 3;
                unsigned a2 = side ? (axis + 2) % 3 : (axis + 1) % 3;
                auto idx = [&](unsigned u, unsigned v) { return (side << axis) | (u << a1) | (v << a2); };
                std::array<unsigned, 4> c{idx(0, 0), idx(1, 0), idx(1, 1), idx(0, 1)};

                // degree-1 bilinear surface, poles u-fastest
                points_vector<T, 3> poles{corner(c[0]), corner(c[1]), corner(c[3]), corner(c[2])};
                std::vector<T> k{0., 0., 1., 1.};
                auto srf = std::make_shared<BSSurface<T, 3>>(poles, k, k, 1, 1);

                std::array<point<T, 2>, 4> uv{{{0., 0.}, {1., 0.}, {1., 1.}, {0., 1.}}};
                Wire<T> w;
                for (unsigned s = 0; s < 4; ++s)
                {
                    unsigned s2 = (s + 1) % 4;
                    Orientation o;
                    auto eid = find_edge(m, b.e, b.v[c[s]], b.v[c[s2]], o);
                    // pcurve parametrized like the edge: from the edge's v1 to its v2
                    auto p_from = o == Orientation::Forward ? uv[s] : uv[s2];
                    auto p_to = o == Orientation::Forward ? uv[s2] : uv[s];
                    auto pc = std::make_shared<BSCurve<T, 2>>(build_segment(p_from, p_to, true));
                    w.coedges.push_back(CoEdge<T>{eid, o, pc});
                }
                w.closed = true;
                auto wid = m.add(std::move(w));
                b.f[fi++] = m.add(Face<T>{srf, {wid}, 1e-7, true});
            }

        Shell sh;
        for (auto fid : b.f)
            sh.faces.push_back(FaceUse{fid, Orientation::Forward});
        sh.closed = true;
        b.shell = m.add(std::move(sh));
        b.solid = m.add(Solid{b.shell, {}});
        return b;
    }
}
