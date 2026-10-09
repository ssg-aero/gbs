// Round trip of the native BREP IGES export (gbs-io/iges_brep.h) through
// OpenCascade: the file written by gbs is read back by IGESControl_Reader,
// every face must be valid, OCCT's own sewing must close the faces into a
// shell, and the enclosed volume must match.
#include <doctest_gtest.hpp>
#include <gbs-io/iges_brep.h>
#include <gbs/bscbuild.h>

#include <BRepBuilderAPI_MakeSolid.hxx>
#include <BRep_Tool.hxx>
#include <BRepBuilderAPI_Sewing.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepGProp.hxx>
#include <GProp_GProps.hxx>
#include <IGESControl_Reader.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Shell.hxx>
#include <TopoDS_Solid.hxx>

#include <cmath>
#include <filesystem>
#include <numbers>

using namespace gbs;
using namespace gbs::brep;
using gbs::operator+;
using gbs::operator*;

namespace
{
    using T = double;

    std::shared_ptr<Surface<T, 3>> quad(point<T, 3> a, point<T, 3> b, point<T, 3> c, point<T, 3> d)
    {
        std::vector<T> k{0., 0., 1., 1.};
        return std::make_shared<BSSurface<T, 3>>(points_vector<T, 3>{a, b, c, d}, k, k, 1, 1);
    }

    std::vector<FaceId> box_faces(Model<T> &m)
    {
        auto corner = [](unsigned i) { return point<T, 3>{T(i & 1), T((i >> 1) & 1), T((i >> 2) & 1)}; };
        std::vector<FaceId> f;
        for (unsigned axis = 0; axis < 3; ++axis)
            for (unsigned side = 0; side < 2; ++side)
            {
                const unsigned a1 = (axis + 1) % 3, a2 = (axis + 2) % 3;
                auto idx = [&](unsigned u, unsigned v) { return (side << axis) | (u << a1) | (v << a2); };
                f.push_back(unwrap(make_face(m, quad(corner(idx(0, 0)), corner(idx(1, 0)), corner(idx(0, 1)), corner(idx(1, 1))))));
            }
        return f;
    }

    // exact rational cylinder of radius 1 and height h
    std::shared_ptr<Surface<T, 3>> nurbs_cylinder(T h)
    {
        auto circle = build_circle<T, 3>(1.);
        points_vector<T, 4> poles;
        for (T z : {T(0), h})
            for (auto p : circle.poles())
            {
                p[2] += p[3] * z;
                poles.push_back(p);
            }
        return std::make_shared<BSSurfaceRational<T, 3>>(poles, circle.knotsFlats(), std::vector<T>{0., 0., h, h}, circle.degree(), 1);
    }

    struct OcctResult
    {
        int faces{};
        bool faces_valid{true};
        bool closed{};
        double volume{};
    };

    OcctResult read_back(const std::string &file, double sew_tol)
    {
        OcctResult r;
        IGESControl_Reader reader;
        REQUIRE(reader.ReadFile(file.c_str()) == IFSelect_RetDone);
        reader.TransferRoots();
        TopoDS_Shape shape = reader.OneShape();
        BRepBuilderAPI_Sewing sewing(sew_tol);
        for (TopExp_Explorer ex(shape, TopAbs_FACE); ex.More(); ex.Next())
        {
            ++r.faces;
            r.faces_valid = r.faces_valid && BRepCheck_Analyzer(ex.Current()).IsValid();
            sewing.Add(ex.Current());
        }
        sewing.Perform();
        TopExp_Explorer shells(sewing.SewedShape(), TopAbs_SHELL);
        REQUIRE(shells.More());
        const auto shell = TopoDS::Shell(shells.Current());
        r.closed = BRep_Tool::IsClosed(shell);
        BRepBuilderAPI_MakeSolid mk(shell);
        GProp_GProps props;
        BRepGProp::VolumeProperties(mk.Solid(), props, 1e-9); // adaptive: the fixed Gauss order is coarse on rational faces
        r.volume = std::abs(props.Mass());
        return r;
    }

    std::string tmp_file(const std::string &name)
    {
        return (std::filesystem::temp_directory_path() / name).string();
    }
}

TEST(tests_brep_iges_occt, box_round_trip)
{
    Model<T> m;
    auto so = unwrap(make_solid(m, unwrap(sew(m, box_faces(m))).shells.front()));
    IgesWriter<T> w;
    w.model().SetUnitsFlag(UNIT_MILLIMETER);
    add_brep(w, m, so, "box");
    const auto file = tmp_file("gbs_brep_box_occt.igs");
    w.write(file);

    auto r = read_back(file, 1e-6);
    ASSERT_EQ(r.faces, 6);
    ASSERT_TRUE(r.faces_valid);
    ASSERT_TRUE(r.closed);
    ASSERT_NEAR(r.volume, 1., 1e-6);
}

TEST(tests_brep_iges_occt, rational_cylinder_round_trip)
{
    // rational surface and rational circles: checks the weighted-pole convention end to end
    Model<T> m;
    const T h = 2.;
    auto lateral = unwrap(make_face(m, nurbs_cylinder(h)));
    auto disk = [&](T z) {
        auto rc = build_circle<T, 3>(1.);
        auto poles = rc.poles();
        for (auto &p : poles)
            p[2] += p[3] * z;
        auto c = std::make_shared<BSCurveRational<T, 3>>(poles, rc.knotsFlats(), rc.degree());
        return unwrap(make_face(m, quad({-2., -2., z}, {2., -2., z}, {-2., 2., z}, {2., 2., z}),
                                unwrap(make_wire(m, std::vector{unwrap(make_edge(m, c))}))));
    };
    auto so = unwrap(make_solid(m, unwrap(sew(m, std::vector{lateral, disk(0.), disk(h)}, {.tol = 1e-5})).shells.front()));

    IgesWriter<T> w;
    w.model().SetUnitsFlag(UNIT_MILLIMETER);
    add_brep(w, m, so, "cyl", IgesExportOptions<T>{.scale = 10.});
    const auto file = tmp_file("gbs_brep_cylinder_occt.igs");
    w.write(file);

    auto r = read_back(file, 1e-3);
    ASSERT_EQ(r.faces, 3);
    ASSERT_TRUE(r.faces_valid);
    ASSERT_TRUE(r.closed);
    const double expected = std::numbers::pi * h * 1000.; // scale 10: lengths x10, volume x1000
    ASSERT_NEAR(r.volume, expected, 1e-6 * expected);
}
