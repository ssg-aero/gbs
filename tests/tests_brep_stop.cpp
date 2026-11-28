// Cooperative stop of the long operations (`StopRequested`, `BuildErrc::Cancelled`): asked
// before, during and after each of them, the model left as it was; the delay of a stop asked
// from another thread measured on heavy cases built by program.
#include <doctest_gtest.hpp>
#include <gbs-brep/blend.h>
#include <gbs-brep/boolean.h>
#include <gbs-brep/equal.h>
#include <gbs-brep/mesh.h>
#include <gbs-brep/pipe.h>
#include <gbs-brep/section.h>
#include <gbs-brep/split.h>
#include <gbs-brep/sweep.h>
#include <gbs-io/step/reader.h>
#include <gbs-io/step/writer.h>
#include <blade_fixture.h>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <functional>
#include <memory>
#include <numbers>
#include <optional>
#include <stdexcept>
#include <ostream>
#include <string>
#include <thread>
#include <vector>

#ifdef GBS_USE_MODULES
    import vecop;
#endif

using namespace gbs;
using namespace gbs::brep;

namespace
{
    using T = double;
    using Clock = std::chrono::steady_clock;
    constexpr T pi = std::numbers::pi_v<T>;
    const point<T, 3> X{1., 0., 0.}, Y{0., 1., 0.}, Z{0., 0., 1.};

    // ---- shapes

    /// A planar face on the polygon p (its plane through p[0], normal along (p1 - p0) x (p2 - p0)).
    FaceId polygon(Model<T> &m, const std::vector<point<T, 3>> &p)
    {
        std::vector<EdgeId> e;
        for (std::size_t i = 0; i < p.size(); ++i)
            e.push_back(unwrap(make_edge(m, p[i], p[(i + 1) % p.size()])));
        auto x = p[1] - p[0];
        x = x / norm(x);
        auto n = cross(p[1] - p[0], p[2] - p[0]);
        n = n / norm(n);
        return unwrap(make_face(m, std::make_shared<Plane<T>>(ax2<T, 3>{{p[0], n, x}}, std::array<T, 4>{-100., 100., -100., 100.}),
                                unwrap(make_wire(m, e))));
    }

    /// The disk of radius r centred at o, normal z.
    FaceId disk(Model<T> &m, T r, const point<T, 3> &o, const point<T, 3> &z = Z, const point<T, 3> &x = X)
    {
        const ax2<T, 3> ax{{o, z, x}};
        const auto e = unwrap(make_edge(m, std::make_shared<EllipseCurve<T>>(r, r, ax)));
        return unwrap(make_face(m, std::make_shared<Plane<T>>(ax, std::array<T, 4>{-2 * r, 2 * r, -2 * r, 2 * r}),
                                unwrap(make_wire(m, std::vector<EdgeId>{e}))));
    }

    SolidId solid_of(Model<T> &m, const BuildResult<SweepResult<T>> &r) { return explore<SolidId>(m, unwrap(r).shape).front(); }

    /// The box [o, o + (a, b, c)].
    SolidId box(Model<T> &m, T a, T b, T c, const point<T, 3> &o = {})
    {
        return solid_of(m, make_prism(m, polygon(m, {o, o + point<T, 3>{a, 0., 0.}, o + point<T, 3>{a, b, 0.}, o + point<T, 3>{0., b, 0.}}),
                                      point<T, 3>{0., 0., c}));
    }

    /// The cylinder of radius r and height h on the disk centred at o.
    SolidId cylinder(Model<T> &m, T r, T h, const point<T, 3> &o = {})
    {
        return solid_of(m, make_prism(m, disk(m, r, o), point<T, 3>{0., 0., h}));
    }

    /// The sphere of radius r centred at o: a half disk turned round the axis through o along z.
    SolidId sphere(Model<T> &m, T r, const point<T, 3> &o)
    {
        const ax2<T, 3> ax{{o, -1. * Y, X}};
        const auto arc = unwrap(make_edge(m, std::make_shared<EllipseCurve<T>>(r, r, ax, -pi / 2, pi / 2)));
        const auto diameter = unwrap(make_edge(m, o + point<T, 3>{0., 0., r}, o - point<T, 3>{0., 0., r}));
        const auto f = unwrap(make_face(m, std::make_shared<Plane<T>>(ax, std::array<T, 4>{-2 * r, 2 * r, -2 * r, 2 * r}),
                                        unwrap(make_wire(m, std::vector<EdgeId>{arc, diameter}))));
        return solid_of(m, make_revolution(m, f, ax1<T, 3>{{o, Z}}));
    }

    /// The edge of the solid whose middle is nearest p.
    EdgeId edge_near(const Model<T> &m, SolidId s, const point<T, 3> &p)
    {
        EdgeId best{};
        T d = std::numeric_limits<T>::max();
        for (auto e : explore<EdgeId>(m, ShapeId{s}))
        {
            const auto &ed = m.edge(e);
            if (ed.degenerate)
                continue;
            const T x = norm(ed.curve->value((ed.u1 + ed.u2) / 2.) - p);
            if (x < d)
                d = x, best = e;
        }
        return best;
    }

    /// The heavy case: a plate of n x n cells pierced by a pin in each cell (cylinders through it) and
    /// spheres half sunk in it between them; `tool` is the compound of the pins and of the spheres.
    struct Plate
    {
        SolidId plate;
        ShapeId tool;
    };
    Plate pierced_plate(Model<T> &m, int n)
    {
        Plate p;
        p.plate = box(m, T(n), T(n), 1.);
        std::vector<ShapeId> tools;
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
            {
                tools.push_back(cylinder(m, 0.2, 2., {i + 0.5, j + 0.5, -0.5}));
                if (i + 1 < n && j + 1 < n)
                    tools.push_back(sphere(m, 0.3, {i + 1., j + 1., 1.}));
            }
        p.tool = unwrap(make_compound(m, tools));
        return p;
    }

    // ---- checks

    /// The model `after` has the same entities alive as `before`, equal (`brep::equal`) once the
    /// dead ones (those an operation added and erased) are compacted away.
    void expect_intact(const Model<T> &before, const Model<T> &after)
    {
        auto a = before, b = after;
        a.compact();
        b.compact();
        const auto d = first_difference(a, b);
        INFO((d ? *d : std::string{}));
        REQUIRE_FALSE(d.has_value());
    }

    template <typename R>
    void expect_cancelled(const R &r)
    {
        REQUIRE_FALSE(r.has_value());
        INFO(r.error().message);
        REQUIRE(r.error().code == BuildErrc::Cancelled);
        REQUIRE(r.error().shapes.empty());
    }

    /**
     * The stop asked of the operation `op(m, stop)` before it starts, at each of its stop points
     * (all of them up to 100, a sample of 50 beyond), and after it ends: the operation fails
     * with `Cancelled`, the model as it was, until the stop is not asked; then it succeeds, as
     * without stop. Returns the number of stop points of a complete run.
     */
    template <typename Op>
    std::size_t expect_stops(Model<T> &m, Op op)
    {
        const auto before = m;
        // before it starts
        expect_cancelled(op(m, [] { return true; }));
        expect_intact(before, m);
        // a complete run, its stop points counted
        std::size_t calls{};
        {
            auto r = op(m, [&] {
                ++calls;
                return false;
            });
            REQUIRE(r.has_value());
        }
        REQUIRE(calls > 0);
        m = before;
        // during: from the k-th stop point on
        std::vector<std::size_t> ks;
        const std::size_t step = calls <= 100 ? 1 : calls / 50;
        for (std::size_t k = 1; k <= calls; k += step)
            ks.push_back(k);
        ks.push_back(calls);
        for (auto k : ks)
        {
            INFO("stop from the stop point " << k << " of " << calls);
            std::size_t c{};
            expect_cancelled(op(m, [&] { return ++c >= k; }));
            expect_intact(before, m);
        }
        // a callback that throws, halfway: a request to stop
        {
            std::size_t c{};
            expect_cancelled(op(m, [&]() -> bool {
                if (++c >= (calls + 1) / 2)
                    throw std::runtime_error("callback failed");
                return false;
            }));
            expect_intact(before, m);
        }
        // after it ends (never asked during the run), then without stop: the operation succeeds
        bool done = false;
        {
            auto r = op(m, [&] { return done; });
            REQUIRE(r.has_value());
            done = true;
        }
        m = before;
        REQUIRE(op(m, {}).has_value());
        return calls;
    }

    /**
     * The operation `op(stop)` run in a thread, the stop asked from this one once it has reached
     * its first stop point and `head` has passed: the delay from the request to the end of the
     * operation (cancelled) in seconds.
     */
    template <typename Op>
    double stop_delay(Op op, std::chrono::milliseconds head = std::chrono::milliseconds{20})
    {
        std::atomic<bool> started{false}, asked{false};
        StopRequested stop = [&] {
            started.store(true, std::memory_order_relaxed);
            return asked.load(std::memory_order_relaxed);
        };
        BuildErrc code{};
        bool failed = false;
        Clock::time_point end;
        std::thread worker{[&] {
            const auto r = op(stop);
            end = Clock::now();
            failed = !r.has_value();
            if (failed)
                code = r.error().code;
        }};
        while (!started.load(std::memory_order_relaxed))
            std::this_thread::yield();
        std::this_thread::sleep_for(head);
        const auto t0 = Clock::now();
        asked.store(true, std::memory_order_relaxed);
        worker.join();
        REQUIRE(failed); // (else the case is too light to measure a stop)
        REQUIRE(code == BuildErrc::Cancelled);
        return std::chrono::duration<double>(end - t0).count();
    }

    /// A complete run of an operation: its time, its stop points, the longest time between two of them
    /// (from its start to the first and from the last to its end included: the delay of a stop at worst).
    struct Run
    {
        double time{}, gap{};
        std::size_t points{};
    };
    std::ostream &operator<<(std::ostream &os, const Run &r)
    {
        return os << r.time << " s, " << r.points << " stop points, the longest time between two " << r.gap << " s";
    }

    template <typename Op>
    Run longest_gap(Op op)
    {
        auto last = Clock::now();
        const auto start = last;
        double gap{};
        std::size_t points{};
        StopRequested stop = [&] {
            ++points;
            const auto now = Clock::now();
            gap = std::max(gap, std::chrono::duration<double>(now - last).count());
            last = now;
            return false;
        };
        REQUIRE(op(stop).has_value());
        const auto end = Clock::now();
        gap = std::max(gap, std::chrono::duration<double>(end - last).count());
        return {std::chrono::duration<double>(end - start).count(), gap, points};
    }

    /// The text of a STEP file of the shapes.
    std::string step_text(const Model<T> &m, const std::vector<ShapeId> &shapes)
    {
        return unwrap(step::write_step(m, std::span<const ShapeId>{shapes})).first;
    }
}

TEST(tests_brep_stop, no_stop_asked)
{
    // an empty stop is never called; a stop never asked changes nothing to the result
    Model<T> m;
    const auto a = box(m, 2., 2., 2.);
    const auto b = cylinder(m, 0.5, 4., {1., 1., -1.});
    auto m1 = m, m2 = m;
    const auto r1 = unwrap(cut(m1, a, b));
    BooleanOptions<T> o;
    o.stop = [] { return false; };
    const auto r2 = unwrap(cut(m2, a, b, o));
    ASSERT_EQ(r1.shape, r2.shape);
    expect_intact(m1, m2);
    ASSERT_EQ(std::string{to_string(BuildErrc::Cancelled)}, "cancelled");
}

TEST(tests_brep_stop, boolean)
{
    Model<T> m;
    const auto a = box(m, 2., 2., 2.);
    const auto b = cylinder(m, 0.5, 4., {1., 1., -1.});
    const auto c = sphere(m, 0.8, {2., 2., 2.});
    for (auto op : {BooleanOp::Fuse, BooleanOp::Cut, BooleanOp::Common})
        for (const auto &tool : {ShapeId{b}, ShapeId{c}})
        {
            const auto n = expect_stops(m, [&](Model<T> &x, StopRequested stop) {
                BooleanOptions<T> o;
                o.stop = std::move(stop);
                return boolean(x, op, a, tool, o);
            });
            MESSAGE("boolean: " << n << " stop points");
        }
}

TEST(tests_brep_stop, section_and_split)
{
    Model<T> m;
    const auto a = box(m, 2., 2., 2.);
    const auto b = sphere(m, 0.8, {2., 2., 2.});
    expect_stops(m, [&](Model<T> &x, StopRequested stop) {
        SectionOptions<T> o;
        o.stop = std::move(stop);
        return section(x, a, b, o);
    });
    const auto plane = std::make_shared<Plane<T>>(ax2<T, 3>{{{0., 0., 1.}, Z, X}}, std::array<T, 4>{-5., 5., -5., 5.});
    expect_stops(m, [&](Model<T> &x, StopRequested stop) {
        SectionOptions<T> o;
        o.stop = std::move(stop);
        return section(x, a, plane, std::optional<std::array<T, 4>>{}, o);
    });
    expect_stops(m, [&](Model<T> &x, StopRequested stop) {
        SectionOptions<T> o;
        o.stop = std::move(stop);
        return split(x, a, b, o);
    });
}

TEST(tests_brep_stop, mesh)
{
    Model<T> m;
    const auto a = box(m, 2., 2., 2.);
    const auto b = sphere(m, 0.8, {2., 2., 2.});
    const auto s = unwrap(fuse(m, a, b)).shape;
    const auto n = expect_stops(m, [&](Model<T> &x, StopRequested stop) {
        MeshOptions<T> o;
        o.stop = std::move(stop);
        return mesh(x, s, o);
    });
    MESSAGE("mesh: " << n << " stop points");
    // a stop inside the refinement of a face: the builder gives Cancelled for that face, and meshes
    // it afterwards
    MeshOptions<T> o;
    bool stop = true;
    o.stop = [&] { return stop; };
    MeshBuilder<T> mb{m, s, o};
    const auto faces = explore<FaceId>(m, s);
    const auto f = *std::ranges::find_if(faces, [&](FaceId x) { return dynamic_cast<const SphericalSurface<T> *>(m.face(x).surface.get()); });
    expect_cancelled(mb.face(f));
    stop = false;
    ASSERT_TRUE(mb.face(f).has_value());
}

TEST(tests_brep_stop, read_step)
{
    // read into a model that already holds shapes: they stay as they were
    Model<T> src;
    const auto a = box(src, 2., 2., 2.);
    const auto b = sphere(src, 0.8, {2., 2., 2.});
    const auto c = cylinder(src, 0.5, 3., {5., 0., 0.});
    const auto text = step_text(src, {unwrap(fuse(src, a, b)).shape, ShapeId{c}});
    Model<T> m;
    box(m, 1., 1., 1., {-3., 0., 0.});
    const auto n = expect_stops(m, [&](Model<T> &x, StopRequested stop) {
        step::StepReadOptions o;
        o.stop = std::move(stop);
        return step::read_step(x, text, o);
    });
    MESSAGE("read_step: " << n << " stop points");
    // strict or not
    expect_stops(m, [&](Model<T> &x, StopRequested stop) {
        step::StepReadOptions o;
        o.strict = true;
        o.stop = std::move(stop);
        return step::read_step(x, text, o);
    });
}

TEST(tests_brep_stop, sweeps)
{
    Model<T> m;
    const auto square = polygon(m, {{0., 0., 0.}, {1., 0., 0.}, {1., 1., 0.}, {0., 1., 0.}});
    expect_stops(m, [&](Model<T> &x, StopRequested stop) {
        SweepOptions<T> o;
        o.stop = std::move(stop);
        return make_prism(x, square, point<T, 3>{0., 0., 2.}, o);
    });
    const auto rect = polygon(m, {{2., 0., 0.}, {3., 0., 0.}, {3., 0., 1.}, {2., 0., 1.}});
    for (const T angle : {pi, 2 * pi})
        expect_stops(m, [&](Model<T> &x, StopRequested stop) {
            SweepOptions<T> o;
            o.stop = std::move(stop);
            return make_revolution(x, rect, ax1<T, 3>{{{0., 0., 0.}, Z}}, angle, o);
        });
    // a disk along an elbow: segment, quarter of circle, segment
    const auto d = disk(m, 0.5, {0., 0., 0.});
    const auto e1 = unwrap(make_edge(m, point<T, 3>{0., 0., 0.}, point<T, 3>{0., 0., 3.}));
    const auto arc = std::make_shared<EllipseCurve<T>>(2., 2., ax2<T, 3>{{{2., 0., 3.}, {0., 1., 0.}, {-1., 0., 0.}}}, 0., pi / 2);
    const auto e2 = unwrap(make_edge(m, arc));
    const auto e3 = unwrap(make_edge(m, arc->value(pi / 2), point<T, 3>{6., 0., 5.}));
    const auto elbow = unwrap(make_wire(m, std::vector<EdgeId>{e1, e2, e3}));
    expect_stops(m, [&](Model<T> &x, StopRequested stop) {
        PipeOptions<T> o;
        o.stop = std::move(stop);
        return make_pipe(x, d, elbow, o);
    });
}

TEST(tests_brep_stop, blends)
{
    Model<T> m;
    const auto s = box(m, 2., 1., 1.5);
    const std::vector<EdgeId> edges{edge_near(m, s, {1., 0., 1.5}), edge_near(m, s, {2., 0., 0.75}), edge_near(m, s, {2., 0.5, 1.5})};
    const auto n = expect_stops(m, [&](Model<T> &x, StopRequested stop) {
        FilletOptions<T> o;
        o.stop = std::move(stop);
        return fillet(x, s, edges, 0.2, o);
    });
    MESSAGE("fillet: " << n << " stop points");
    expect_stops(m, [&](Model<T> &x, StopRequested stop) {
        FilletOptions<T> o;
        o.stop = std::move(stop);
        return chamfer(x, s, edges, 0.1, o);
    });
}

TEST(tests_brep_stop, delays)
{
    // heavy cases: a plate of 6 x 6 cells, its pins and spheres (a boolean of 62 solids, a mesh of its
    // result, the reading of a file of the operands and of five copies of the result); a stop asked from another
    // thread ends the operation in less than a second, the model as it was, and the operation goes
    // through afterwards. The longest time between two stop points is the delay of a stop at worst.
    Model<T> m;
    const auto p = pierced_plate(m, 6);
    const auto before = m;
    BooleanOptions<T> bo;
    const double db = stop_delay([&](StopRequested stop) {
        bo.stop = std::move(stop);
        return cut(m, p.plate, p.tool, bo);
    });
    expect_intact(before, m);
    std::optional<BooleanResult<T>> r;
    const auto run_b = longest_gap([&](StopRequested stop) {
        bo.stop = std::move(stop);
        auto x = cut(m, p.plate, p.tool, bo);
        r = x ? std::optional{*x} : std::nullopt;
        return x;
    });
    MESSAGE("boolean: " << run_b << "; stopped in " << db << " s");
    CHECK(db < 1.);
    CHECK(run_b.gap < 1.);

    MeshOptions<T> mo;
    mo.deflection = 1e-4;
    const auto run_m = longest_gap([&](StopRequested stop) {
        mo.stop = std::move(stop);
        return mesh(m, r->shape, mo);
    });
    const double dm = stop_delay([&](StopRequested stop) {
        mo.stop = std::move(stop);
        return mesh(m, r->shape, mo);
    });
    MESSAGE("mesh: " << run_m << "; stopped in " << dm << " s");
    CHECK(dm < 1.);
    CHECK(run_m.gap < 1.);

    std::vector<ShapeId> shapes{ShapeId{p.plate}, p.tool, r->shape};
    for (int k = 1; k <= 4; ++k) // (copies of the result, for a reading long enough to be stopped)
        shapes.push_back(copy_shape(m, r->shape, Placement<T>{.t = {0., 0., 2. * k}}));
    const auto text = step_text(m, shapes);
    Model<T> target;
    box(target, 1., 1., 1., {-3., 0., 0.});
    const auto target_before = target;
    step::StepReadOptions so;
    const double ds = stop_delay([&](StopRequested stop) {
        so.stop = std::move(stop);
        return step::read_step(target, text, so);
    });
    expect_intact(target_before, target);
    const auto run_s = longest_gap([&](StopRequested stop) {
        so.stop = std::move(stop);
        return step::read_step(target, text, so);
    });
    MESSAGE("read_step: " << text.size() / 1024 << " kB, " << run_s << "; stopped in " << ds << " s");
    CHECK(ds < 1.);
    CHECK(run_s.gap < 1.);
}

TEST(tests_brep_stop, delays_bspline)
{
    // a blade skinned through 24 twisted sections (one B-spline face of many patches) cut by a copy of
    // it twisted and moved a little: the faces lie close over their whole extent, and the search of the
    // partial coincidences (the first step of the boolean, `split_partial_coincidences`) examines
    // every patch of each against the other
    Model<T> m;
    std::vector<T> zs, twists, twists2;
    for (int i = 0; i < 24; ++i)
    {
        zs.push_back(2. * i / 23.);
        twists.push_back(-0.4 + 0.8 * i / 23. + 0.05 * std::sin(3. * i));
        twists2.push_back(twists.back() + 0.02);
    }
    const auto a = ShapeId{blade_fixture::blade(m, zs, twists)};
    const auto b = copy_shape(m, ShapeId{blade_fixture::blade(m, zs, twists2)}, Placement<T>{.t = {0.01, 0.005, 0.3}});
    const auto before = m;
    BooleanOptions<T> bo;
    const double db = stop_delay([&](StopRequested stop) {
        bo.stop = std::move(stop);
        return cut(m, a, b, bo);
    });
    expect_intact(before, m);
    const auto run_b = longest_gap([&](StopRequested stop) {
        bo.stop = std::move(stop);
        return cut(m, a, b, bo);
    });
    MESSAGE("boolean of B-spline faces: " << run_b << "; stopped in " << db << " s");
    CHECK(db < 1.);
    CHECK(run_b.gap < 1.);
}
