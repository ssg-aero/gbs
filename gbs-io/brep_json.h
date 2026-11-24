#pragma once
/**
 * @file brep_json.h
 * @brief A BREP model (`gbs::brep::Model<T>`) written to JSON and read back exactly.
 *
 * Design: docs/sources/design/io_brep_json.md.
 *
 * What is read back is the model that was written, not a copy of its shapes:
 * - every table (vertices, edges, wires, faces, shells, solids, compounds) at the same size, the
 *   dead entities included (`null` in the file, dead again once read): every id stays valid, an
 *   id held outside the model (an application's document) designates the same entity. Calling
 *   `Model::compact()` before writing drops them, its `IdRemap` translating the ids held outside;
 * - the ids stored in the entities as they are, an invalid id (npos) and an id of a dead entity
 *   included (`erase` does not cascade); one beyond its table is refused when writing;
 * - orientations, tolerances, flags, names, external ids, instances of parts, the next part key,
 *   the unit;
 * - the geometry, written once per object and referred to by its index in the file: a curve or a
 *   surface shared by several entities (or inside another one: the basis of an offset, the
 *   surface of a curve on surface) is shared again once read, an object used once is not;
 * - every number bit for bit (the shortest decimal form that reads back the same double, read
 *   with RapidJSON's kParseFullPrecisionFlag; -0, the subnormals and the infinities included; a
 *   NaN comes back a NaN, without its sign and payload).
 *
 * The geometry written is that of gbs's own classes, recognized by their exact dynamic type:
 * B-spline and NURBS curves and surfaces (2D and 3D), the primitives (line, ellipse, plane,
 * cylinder, cone, sphere, torus), curves trimmed, composite, reversed, reparametrized, extended,
 * on a surface, offset (by a B-spline function), moved (`brep::MovedCurve`), surfaces of
 * revolution, offset and moved. Another class (a curve of an application, the old `Line` and
 * `Circle` whose constructors normalize what they are given) cannot be read back exactly: the
 * writer refuses it (`BrepWriteError`) instead of writing something else.
 *
 * Reading refuses, with `BrepReadError` (its code, the JSON pointer of the value at fault and the
 * reason): a text that is not JSON or is cut short, a JSON text that is not a model of gbs, a
 * version newer than this code, and anything that a model could not hold: a key unknown or
 * repeated, a value missing or of the wrong type, an index beyond its table, a geometry of the
 * wrong kind or referring to one after it, the knots of a B-spline out of order or of the wrong
 * number. Nothing in the text is trusted before it is checked: the text is parsed without
 * recursion (a deep nesting does not exhaust the stack), sizes come from the arrays read, and the
 * constructors of gbs check the rest. What is not checked is what the model itself does not
 * check: a model read is one the `Model` API could have built, not a valid solid (`check` says).
 */
#include <gbs-io/json.h>
#include <gbs-brep/model.h>
#include <gbs-brep/edit.h>
#include <gbs-brep/seamless.h>
#include <gbs/curves>
#include <gbs/surfaces>
#include <gbs/primitives.h>
#include <gbs/offsets.h>

#include <rapidjson/document.h>
#include <rapidjson/prettywriter.h>
#include <rapidjson/stringbuffer.h>
#include <rapidjson/writer.h>

#include <array>
#include <bit>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <format>
#include <iterator>
#include <istream>
#include <memory>
#include <optional>
#include <ostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <typeinfo>
#include <unordered_map>
#include <utility>
#include <variant>
#include <vector>

namespace gbs::io
{
    /// The value of "format" at the top of a file.
    inline constexpr std::string_view brep_json_format = "gbs-brep";
    /// The version written, and the newest one read.
    inline constexpr std::uint32_t brep_json_version = 1;

    /// The reals a model is written in exactly: those of a JSON number (a double).
    template <typename T>
    concept brep_json_real = std::same_as<T, double> || std::same_as<T, float>;

    /// Why a text is not a model that can be read.
    enum class BrepReadErrc : std::uint8_t
    {
        Syntax,       ///< not a JSON text, or one cut short
        NotBrep,      ///< a JSON text, but not a BREP model of gbs ("format")
        NewerVersion, ///< a model written by a newer version of the format
        Invalid,      ///< a value missing, unknown, of the wrong type, out of its range, or inconsistent
    };

    [[nodiscard]] inline constexpr const char *to_string(BrepReadErrc c) noexcept
    {
        switch (c)
        {
        case BrepReadErrc::Syntax: return "syntax";
        case BrepReadErrc::NotBrep: return "not a BREP model";
        case BrepReadErrc::NewerVersion: return "newer version";
        case BrepReadErrc::Invalid: return "invalid";
        }
        return "invalid";
    }

    /**
     * @brief A text that is not a model that can be read: the code, the JSON pointer of the value
     * at fault (`where`, as "/edges/12/u/0"; empty for the whole text) and the reason. what() is
     * "source: where: reason"; where and reason are views of it, so that copying the error never
     * throws.
     */
    class BrepReadError : public TextError
    {
        BrepReadErrc code_;
        std::size_t where_at_{}, where_size_{};

        static auto message(std::string_view source, std::string_view where, std::string_view reason) -> std::string
        {
            return where.empty() ? std::format("{}: {}", source, reason) : std::format("{}: {}: {}", source, where, reason);
        }

    public:
        BrepReadError(std::string_view source, BrepReadErrc code, std::string_view where, std::string_view reason)
            : TextError{message(source, where, reason), reason}, code_{code}, where_at_{where.empty() ? 0 : source.size() + 2}, where_size_{where.size()}
        {
        }
        [[nodiscard]] auto code() const noexcept -> BrepReadErrc { return code_; }
        /// The JSON pointer of the value at fault, empty for the whole text.
        [[nodiscard]] auto where() const noexcept -> std::string_view { return std::string_view{what()}.substr(where_at_, where_size_); }
    };

    /// A model that cannot be written exactly: a geometry of a class the format does not hold, an
    /// id beyond its table. Nothing is written.
    class BrepWriteError : public BRepError
    {
    public:
        explicit BrepWriteError(const std::string &message) : BRepError{"write_brep: " + message} {}
    };

    struct BrepWriteOptions
    {
        bool pretty{false}; ///< indented, an array of numbers on one line; compact otherwise
    };

    namespace brep_json_detail
    {
        using rapidjson::SizeType;
        using rapidjson::Value;

        inline constexpr std::array<const char *, 7> table_keys{"vertices", "edges", "wires", "faces", "shells", "solids", "compounds"};
        inline constexpr std::array<const char *, 7> type_names{"vertex", "edge", "wire", "face", "shell", "solid", "compound"};

        template <typename T>
        constexpr const char *real_name = std::same_as<T, double> ? "double" : "float";

        /// Same bits (0 and -0 differ, a NaN equals itself).
        template <std::floating_point T>
        [[nodiscard]] bool same_bits(T a, T b) noexcept
        {
            return std::memcmp(&a, &b, sizeof(T)) == 0;
        }
        template <std::floating_point T, std::size_t N>
        [[nodiscard]] bool same_bits(const std::array<T, N> &a, const std::array<T, N> &b) noexcept
        {
            return std::memcmp(a.data(), b.data(), sizeof(T) * N) == 0;
        }

        /// What gbs's B-spline constructors accept (`check_curve`), the knots finite: the condition
        /// of a B-spline written and read.
        template <typename T>
        [[nodiscard]] bool valid_knots(const std::vector<T> &k, std::size_t poles, std::size_t degree) noexcept
        {
            if (degree >= poles || k.size() != poles + degree + 1)
                return false;
            for (std::size_t i = 0; i < k.size(); ++i)
                if (!std::isfinite(k[i]) || (i > 0 && k[i] < k[i - 1]))
                    return false;
            return true;
        }

        /// U bounds of a B-spline surface that differ from its knots, as `trimU` leaves them when
        /// they are within knot_eps of the ends (it then changes neither knots nor poles).
        template <typename T>
        [[nodiscard]] bool trim_u_noop(const std::vector<T> &k, T u1, T u2) noexcept
        {
            return std::fmin(u1, u2) - k.front() < knot_eps<T> && k.back() - std::fmax(u1, u2) < knot_eps<T>;
        }

        /// The geometry classes of a pool entry.
        enum class Kind : std::uint8_t
        {
            Curve2,
            Curve3,
            Surface2,
            Surface3
        };

        template <typename T, std::size_t D>
        consteval Kind curve_kind() { return D == 2 ? Kind::Curve2 : Kind::Curve3; }
        template <typename T, std::size_t D>
        consteval Kind surface_kind() { return D == 2 ? Kind::Surface2 : Kind::Surface3; }

        inline const char *kind_name(Kind k) noexcept
        {
            switch (k)
            {
            case Kind::Curve2: return "a 2D curve";
            case Kind::Curve3: return "a 3D curve";
            case Kind::Surface2: return "a 2D surface";
            case Kind::Surface3: return "a 3D surface";
            }
            return "a geometry";
        }

        /// `f(static_cast<const X &>(b))` if the dynamic type of b is exactly X.
        template <typename X, typename B, typename F>
        bool as(const B &b, F &f)
        {
            if (typeid(b) != typeid(X))
                return false;
            f(static_cast<const X &>(b));
            return true;
        }

        /// f applied to the curve as its exact class; false for a class the format does not hold.
        template <typename T, std::size_t D, typename F>
        bool visit_curve(const Curve<T, D> &c, F &&f)
        {
            if (as<BSCurve<T, D>>(c, f) || as<BSCurveRational<T, D>>(c, f) || as<CurveTrimmed<T, D>>(c, f) || as<CurveComposite<T, D>>(c, f) ||
                as<CurveReversed<T, D>>(c, f) || as<CurveReparametrized<T, D>>(c, f) || as<CurveExtended<T, D>>(c, f) ||
                as<CurveOnSurface<T, D>>(c, f))
                return true;
            if constexpr (D == 3)
                return as<LineCurve<T>>(c, f) || as<EllipseCurve<T>>(c, f) || as<brep::MovedCurve<T>>(c, f) || as<CurveOffset3D<T, BSCfunction<T>>>(c, f);
            else
                return as<CurveOffset2D<T, BSCfunction<T>>>(c, f);
        }

        /// f applied to the surface as its exact class; false for a class the format does not hold.
        template <typename T, std::size_t D, typename F>
        bool visit_surface(const Surface<T, D> &s, F &&f)
        {
            if (as<BSSurface<T, D>>(s, f) || as<BSSurfaceRational<T, D>>(s, f))
                return true;
            if constexpr (D == 3)
                return as<Plane<T>>(s, f) || as<CylindricalSurface<T>>(s, f) || as<ConicalSurface<T>>(s, f) || as<SphericalSurface<T>>(s, f) ||
                       as<ToroidalSurface<T>>(s, f) || as<SurfaceOfRevolution<T>>(s, f) || as<OffsetSurface<T>>(s, f) ||
                       as<brep::MovedSurface<T>>(s, f) || as<brep::SeamlessSurface<T>>(s, f);
            else
                return false;
        }

        // =====================================================================
        // Writing
        // =====================================================================

        template <typename T, typename W>
        class Writer
        {
            W &w_;
            const brep::Model<T> &m_;

            struct Node
            {
                Kind kind;
                const void *geom; ///< const Curve<T, 2> *, const Curve<T, 3> *, const Surface<T, 2> * or const Surface<T, 3> *
            };
            std::unordered_map<const void *, std::uint32_t> index_; ///< most derived address -> index in the pool
            std::vector<Node> nodes_;                               ///< the pool, each geometry after those it refers to
            std::string ctx_;                                       ///< what is being written, for an error

            [[noreturn]] void fail(std::string_view why) const { throw BrepWriteError{ctx_.empty() ? std::string{why} : std::format("{}: {}", ctx_, why)}; }

            static std::string class_name(const auto &x) { return typeid(x).name(); }

            // ---- pass 1: the pool, children first --------------------------------------------

            template <std::size_t D>
            auto add(const Curve<T, D> *c) -> std::uint32_t
            {
                if (!c)
                    fail("null curve inside a geometry");
                const void *key = dynamic_cast<const void *>(c);
                if (auto it = index_.find(key); it != index_.end())
                    return it->second;
                if (!visit_curve<T, D>(*c, [this](const auto &x) { this->children(x); }))
                    fail(std::format("curve of class {} cannot be written: it is not one of the classes of gbs the format holds", class_name(*c)));
                return push(key, curve_kind<T, D>(), c);
            }

            template <std::size_t D>
            auto add(const Surface<T, D> *s) -> std::uint32_t
            {
                if (!s)
                    fail("null surface inside a geometry");
                const void *key = dynamic_cast<const void *>(s);
                if (auto it = index_.find(key); it != index_.end())
                    return it->second;
                if (!visit_surface<T, D>(*s, [this](const auto &x) { this->children(x); }))
                    fail(std::format("surface of class {} cannot be written: it is not one of the classes of gbs the format holds", class_name(*s)));
                return push(key, surface_kind<T, D>(), s);
            }

            auto push(const void *key, Kind kind, const void *geom) -> std::uint32_t
            {
                const auto i = static_cast<std::uint32_t>(nodes_.size());
                index_.emplace(key, i);
                nodes_.push_back({kind, geom});
                return i;
            }

            template <std::size_t D, bool R>
            void children(const BSCurveGeneral<T, D, R> &c)
            {
                if (!valid_knots(c.knotsFlats(), c.poles().size(), c.degree()))
                    fail("B-spline curve with knots out of order, not finite or of the wrong number");
            }
            template <std::size_t D, bool R>
            void children(const BSSurfaceGeneral<T, D, R> &s)
            {
                if (!valid_knots(s.knotsFlatsU(), s.nPolesU(), s.degreeU()) || !valid_knots(s.knotsFlatsV(), s.nPolesV(), s.degreeV()) ||
                    s.nPolesU() * s.nPolesV() != s.poles().size())
                    fail("B-spline surface with knots out of order, not finite or of the wrong number");
                const auto b = s.bounds();
                const auto &ku = s.knotsFlatsU();
                const auto &kv = s.knotsFlatsV();
                // bounds that the reader can give back: U by trimU (no change of the knots), V not at all
                if (!same_bits(b[2], kv.front()) || !same_bits(b[3], kv.back()))
                    fail("B-spline surface whose V bounds are not its end knots: no constructor gives it back");
                if ((!same_bits(b[0], ku.front()) || !same_bits(b[1], ku.back())) && !trim_u_noop(ku, b[0], b[1]))
                    fail("B-spline surface whose U bounds are more than knot_eps inside its end knots: no constructor gives it back");
            }
            void children(const BSCfunction<T> &f) { children(f.basisCurve()); }
            void children(const LineCurve<T> &) {}
            void children(const EllipseCurve<T> &) {}
            void children(const Plane<T> &) {}
            void children(const CylindricalSurface<T> &) {}
            void children(const ConicalSurface<T> &) {}
            void children(const SphericalSurface<T> &) {}
            void children(const ToroidalSurface<T> &) {}
            template <std::size_t D>
            void children(const CurveTrimmed<T, D> &c) { add<D>(c.basisCurve().get()); }
            template <std::size_t D>
            void children(const CurveComposite<T, D> &c)
            {
                if (c.curves().empty())
                    fail("composite curve of no curve");
                for (const auto &x : c.curves())
                    add<D>(x.get());
            }
            template <std::size_t D>
            void children(const CurveReversed<T, D> &c) { add<D>(c.basisCurve().get()); }
            template <std::size_t D>
            void children(const CurveReparametrized<T, D> &c)
            {
                add<D>(c.basisCurve().get());
                children(c.paramFunction());
            }
            template <std::size_t D>
            void children(const CurveExtended<T, D> &c) { add<D>(c.basisCurve().get()); }
            template <std::size_t D>
            void children(const CurveOnSurface<T, D> &c)
            {
                add<2>(c.p_basisCurve().get());
                add<D>(c.p_basisSurface().get());
            }
            void children(const CurveOffset2D<T, BSCfunction<T>> &c)
            {
                add<2>(&c.basisCurve());
                children(c.offset());
            }
            void children(const CurveOffset3D<T, BSCfunction<T>> &c)
            {
                add<3>(&c.basisCurve());
                children(c.offset());
                const auto d = c.direction();
                if (!same_bits(normalized(d), d)) // the constructor normalizes: d must come back as it is
                    fail("offset curve whose direction is not given back by its normalization");
            }
            void children(const brep::MovedCurve<T> &c) { add<3>(c.basis().get()); }
            void children(const SurfaceOfRevolution<T> &s) { add<2>(s.basisCurve().get()); }
            void children(const OffsetSurface<T> &s) { add<3>(s.basis().get()); }
            void children(const brep::MovedSurface<T> &s) { add<3>(s.basis().get()); }
            void children(const brep::SeamlessSurface<T> &s) { add<3>(s.basis().get()); }

            // ---- pass 2: writing -------------------------------------------------------------

            void key(const char *k) { w_.Key(k); }
            void real(T x) { w_.Double(static_cast<double>(x)); }
            void uint(std::uint64_t x) { w_.Uint64(x); }
            void boolean(bool b) { w_.Bool(b); }
            template <std::size_t N>
            void reals(const std::array<T, N> &a)
            {
                w_.StartArray();
                for (T x : a)
                    real(x);
                w_.EndArray();
            }
            void reals(const std::vector<T> &v)
            {
                w_.StartArray();
                for (T x : v)
                    real(x);
                w_.EndArray();
            }
            template <std::size_t N>
            void points(const std::vector<std::array<T, N>> &v)
            {
                w_.StartArray();
                for (const auto &p : v)
                    reals(p);
                w_.EndArray();
            }
            void frame(const Frame3<T> &f)
            {
                w_.StartArray();
                reals(f.o);
                reals(f.x);
                reals(f.y);
                reals(f.z);
                w_.EndArray();
            }
            void placement(const brep::Placement<T> &p)
            {
                w_.StartObject();
                key("r");
                w_.StartArray();
                for (const auto &row : p.r)
                    reals(row);
                w_.EndArray();
                key("t");
                reals(p.t);
                w_.EndObject();
            }
            /// a reference to a geometry of the pool (keyed by its most derived address)
            template <typename G>
            void ref(const G *g) { w_.Uint(index_.at(dynamic_cast<const void *>(g))); }

            void head(const char *type, std::size_t dim)
            {
                key("type");
                w_.String(type);
                key("dim");
                w_.Uint(static_cast<unsigned>(dim));
            }

            template <std::size_t N>
            void bspline(const BSCurveGeneral<T, N, false> &c) // a function: inline
            {
                w_.StartObject();
                bspline_fields(c);
                w_.EndObject();
            }
            template <std::size_t D, bool R>
            void bspline_fields(const BSCurveGeneral<T, D, R> &c)
            {
                key("degree");
                uint(c.degree());
                key("knots");
                reals(c.knotsFlats());
                key("poles");
                points(c.poles()); // homogeneous (w x, w y, w z, w) for a NURBS, as stored
                const auto b = c.bounds();
                if (!same_bits(b[0], c.knotsFlats().front()) || !same_bits(b[1], c.knotsFlats().back()))
                {
                    key("bounds");
                    reals(b);
                }
            }
            void function(const BSCfunction<T> &f) { bspline(f.basisCurve()); }

            template <std::size_t D>
            void write(const BSCurve<T, D> &c) { head("bspline_curve", D), bspline_fields(c); }
            template <std::size_t D>
            void write(const BSCurveRational<T, D> &c) { head("nurbs_curve", D), bspline_fields(c); }
            template <std::size_t D, bool R>
            void bspline_surface(const BSSurfaceGeneral<T, D, R> &s)
            {
                head(R ? "nurbs_surface" : "bspline_surface", D);
                key("degree");
                w_.StartArray();
                uint(s.degreeU());
                uint(s.degreeV());
                w_.EndArray();
                key("knots");
                w_.StartArray();
                reals(s.knotsFlatsU());
                reals(s.knotsFlatsV());
                w_.EndArray();
                key("poles");
                points(s.poles()); // u first: pole (i, j) at i + j nu
                const auto b = s.bounds();
                if (!same_bits(b[0], s.knotsFlatsU().front()) || !same_bits(b[1], s.knotsFlatsU().back()))
                {
                    key("bounds");
                    reals(b);
                }
            }
            template <std::size_t D>
            void write(const BSSurface<T, D> &s) { bspline_surface(s); }
            template <std::size_t D>
            void write(const BSSurfaceRational<T, D> &s) { bspline_surface(s); }
            void write(const LineCurve<T> &c)
            {
                head("line", 3);
                key("origin");
                reals(c.origin());
                key("vector");
                reals(c.vector());
                key("range");
                reals(c.bounds());
            }
            void write(const EllipseCurve<T> &c)
            {
                head("ellipse", 3);
                key("radii");
                reals(std::array<T, 2>{c.radius1(), c.radius2()});
                key("frame");
                frame(c.frame());
                key("range");
                reals(c.bounds());
            }
            void elementary(const char *type, const ElementarySurface<T> &s, std::initializer_list<std::pair<const char *, T>> params)
            {
                head(type, 3);
                for (const auto &[k, x] : params)
                {
                    key(k);
                    real(x);
                }
                key("frame");
                frame(s.frame());
                key("box");
                reals(s.bounds());
            }
            void write(const Plane<T> &s) { elementary("plane", s, {}); }
            void write(const CylindricalSurface<T> &s) { elementary("cylinder", s, {{"radius", s.radius()}}); }
            void write(const ConicalSurface<T> &s) { elementary("cone", s, {{"radius", s.radius()}, {"semi_angle", s.semi_angle()}}); }
            void write(const SphericalSurface<T> &s) { elementary("sphere", s, {{"radius", s.radius()}}); }
            void write(const ToroidalSurface<T> &s) { elementary("torus", s, {{"major_radius", s.major_radius()}, {"minor_radius", s.minor_radius()}}); }
            template <std::size_t D>
            void write(const CurveTrimmed<T, D> &c)
            {
                head("trimmed_curve", D);
                key("basis");
                ref(c.basisCurve().get());
                key("range");
                reals(c.bounds());
            }
            template <std::size_t D>
            void write(const CurveComposite<T, D> &c)
            {
                head("composite_curve", D);
                key("curves");
                w_.StartArray();
                for (const auto &x : c.curves())
                    ref(x.get());
                w_.EndArray();
            }
            template <std::size_t D>
            void write(const CurveReversed<T, D> &c)
            {
                head("reversed_curve", D);
                key("basis");
                ref(c.basisCurve().get());
            }
            template <std::size_t D>
            void write(const CurveReparametrized<T, D> &c)
            {
                head("reparametrized_curve", D);
                key("basis");
                ref(c.basisCurve().get());
                key("function");
                function(c.paramFunction());
            }
            template <std::size_t D>
            void write(const CurveExtended<T, D> &c)
            {
                head("extended_curve", D);
                key("basis");
                ref(c.basisCurve().get());
                key("clamped");
                boolean(c.isClamped());
            }
            template <std::size_t D>
            void write(const CurveOnSurface<T, D> &c)
            {
                head("curve_on_surface", D);
                key("pcurve");
                ref(c.p_basisCurve().get());
                key("surface");
                ref(c.p_basisSurface().get());
            }
            void write(const CurveOffset2D<T, BSCfunction<T>> &c)
            {
                head("offset_curve", 2);
                key("basis");
                ref(&c.basisCurve());
                key("function");
                function(c.offset());
            }
            void write(const CurveOffset3D<T, BSCfunction<T>> &c)
            {
                head("offset_curve", 3);
                key("basis");
                ref(&c.basisCurve());
                key("function");
                function(c.offset());
                key("direction");
                reals(c.direction());
            }
            void write(const brep::MovedCurve<T> &c)
            {
                head("moved_curve", 3);
                key("basis");
                ref(c.basis().get());
                key("placement");
                placement(c.placement());
            }
            void write(const SurfaceOfRevolution<T> &s)
            {
                head("revolution_surface", 3);
                key("curve");
                ref(s.basisCurve().get());
                key("axis");
                w_.StartArray();
                for (const auto &p : s.axis2())
                    reals(p);
                w_.EndArray();
                const auto b = s.bounds();
                key("angles");
                reals(std::array<T, 2>{b[2], b[3]});
            }
            void write(const OffsetSurface<T> &s)
            {
                head("offset_surface", 3);
                key("basis");
                ref(s.basis().get());
                key("distance");
                real(s.distance());
            }
            void write(const brep::MovedSurface<T> &s)
            {
                head("moved_surface", 3);
                key("basis");
                ref(s.basis().get());
                key("placement");
                placement(s.placement());
            }

            void write(const brep::SeamlessSurface<T> &s)
            {
                head("seamless_surface", 3);
                key("basis");
                ref(s.basis().get());
                key("periods");
                reals(s.periods());
            }

            void write_node(const Node &n)
            {
                w_.StartObject();
                const auto f = [this](const auto &x) { this->write(x); };
                switch (n.kind)
                {
                case Kind::Curve2: visit_curve<T, 2>(*static_cast<const Curve<T, 2> *>(n.geom), f); break;
                case Kind::Curve3: visit_curve<T, 3>(*static_cast<const Curve<T, 3> *>(n.geom), f); break;
                case Kind::Surface2: visit_surface<T, 2>(*static_cast<const Surface<T, 2> *>(n.geom), f); break;
                case Kind::Surface3: visit_surface<T, 3>(*static_cast<const Surface<T, 3> *>(n.geom), f); break;
                }
                w_.EndObject();
            }

            // ---- topology --------------------------------------------------------------------

            /// An id stored in an entity: null if invalid, its index otherwise, below the size of its table.
            template <brep::ShapeType K>
            void id(brep::Id<K> h)
            {
                if (!h.valid())
                {
                    w_.Null();
                    return;
                }
                if (h.index >= m_.template capacity<K>())
                    fail(std::format("{} id {} beyond its table ({} {}s)", type_names[std::to_underlying(K)], h.index, m_.template capacity<K>(),
                                     type_names[std::to_underlying(K)]));
                w_.Uint(h.index);
            }
            void id(const brep::ShapeId &s)
            {
                w_.StartArray();
                w_.String(type_names[std::to_underlying(brep::shape_type(s))]);
                std::visit([this](auto h) { this->id(h); }, s);
                w_.EndArray();
            }
            void orient(brep::Orientation o) { w_.String(o == brep::Orientation::Forward ? "F" : "R"); }

            /// name, external id, instance of a live entity, if it has them
            void attributes(const brep::ShapeId &s)
            {
                if (auto n = m_.name(s); !n.empty())
                {
                    key("name");
                    w_.String(n.data(), static_cast<SizeType>(n.size()));
                }
                if (auto x = m_.externalId(s))
                {
                    key("external_id");
                    w_.Int64(*x);
                }
                if (auto i = m_.instance(s))
                {
                    key("instance");
                    w_.StartObject();
                    key("part");
                    uint(i->part);
                    key("placement");
                    placement(i->place);
                    w_.EndObject();
                }
            }

            void entity(const brep::Vertex<T> &v)
            {
                key("p");
                reals(v.pnt);
                key("tol");
                real(v.tol);
            }
            void entity(const brep::Edge<T> &e)
            {
                key("curve");
                if (e.curve)
                    ref(e.curve.get());
                else
                    w_.Null();
                key("u");
                reals(std::array<T, 2>{e.u1, e.u2});
                key("v");
                w_.StartArray();
                id(e.v1);
                id(e.v2);
                w_.EndArray();
                key("tol");
                real(e.tol);
                key("degenerate");
                boolean(e.degenerate);
                key("same_parameter");
                boolean(e.same_parameter);
            }
            void entity(const brep::Wire<T> &w)
            {
                key("coedges");
                w_.StartArray();
                for (const auto &ce : w.coedges)
                {
                    w_.StartArray();
                    id(ce.edge);
                    orient(ce.orient);
                    if (ce.pcurve)
                        ref(ce.pcurve.get());
                    else
                        w_.Null();
                    w_.EndArray();
                }
                w_.EndArray();
                key("closed");
                boolean(w.closed);
            }
            void entity(const brep::Face<T> &f)
            {
                key("surface");
                if (f.surface)
                    ref(f.surface.get());
                else
                    w_.Null();
                key("wires");
                w_.StartArray();
                for (auto h : f.wires)
                    id(h);
                w_.EndArray();
                key("tol");
                real(f.tol);
                key("natural_bounds");
                boolean(f.natural_bounds);
            }
            void entity(const brep::Shell &s)
            {
                key("faces");
                w_.StartArray();
                for (const auto &fu : s.faces)
                {
                    w_.StartArray();
                    id(fu.face);
                    orient(fu.orient);
                    w_.EndArray();
                }
                w_.EndArray();
                key("closed");
                boolean(s.closed);
            }
            void entity(const brep::Solid &s)
            {
                key("outer");
                id(s.outer);
                key("voids");
                w_.StartArray();
                for (auto h : s.voids)
                    id(h);
                w_.EndArray();
            }
            void entity(const brep::Compound &c)
            {
                key("shapes");
                w_.StartArray();
                for (const auto &s : c.shapes)
                    id(s);
                w_.EndArray();
            }

            template <brep::ShapeType K>
            void table()
            {
                key(table_keys[std::to_underlying(K)]);
                w_.StartArray();
                const auto n = m_.template capacity<K>();
                for (std::uint32_t i = 0; i < n; ++i)
                {
                    const brep::Id<K> h{i};
                    if (!m_.alive(h))
                    {
                        w_.Null();
                        continue;
                    }
                    ctx_ = std::format("{} {}", type_names[std::to_underlying(K)], i);
                    w_.StartObject();
                    if constexpr (K == brep::ShapeType::Vertex) entity(m_.vertex(h));
                    else if constexpr (K == brep::ShapeType::Edge) entity(m_.edge(h));
                    else if constexpr (K == brep::ShapeType::Wire) entity(m_.wire(h));
                    else if constexpr (K == brep::ShapeType::Face) entity(m_.face(h));
                    else if constexpr (K == brep::ShapeType::Shell) entity(m_.shell(h));
                    else if constexpr (K == brep::ShapeType::Solid) entity(m_.solid(h));
                    else entity(m_.compound(h));
                    attributes(h);
                    w_.EndObject();
                }
                w_.EndArray();
            }

            /// pass 1: the geometry of the live entities, in the order of the tables
            void collect()
            {
                for (auto h : m_.template ids<brep::EdgeId>())
                {
                    ctx_ = std::format("edge {}", h.index);
                    if (const auto &c = m_.edge(h).curve)
                        add<3>(c.get());
                }
                for (auto h : m_.template ids<brep::WireId>())
                {
                    const auto &w = m_.wire(h);
                    for (std::size_t k = 0; k < w.coedges.size(); ++k)
                    {
                        ctx_ = std::format("wire {}, co-edge {}", h.index, k);
                        if (const auto &c = w.coedges[k].pcurve)
                            add<2>(c.get());
                    }
                }
                for (auto h : m_.template ids<brep::FaceId>())
                {
                    ctx_ = std::format("face {}", h.index);
                    if (const auto &s = m_.face(h).surface)
                        add<3>(s.get());
                }
                ctx_.clear();
            }

        public:
            Writer(W &w, const brep::Model<T> &m) : w_{w}, m_{m} {}

            void run()
            {
                collect();
                w_.StartObject();
                key("format");
                w_.String(brep_json_format.data(), static_cast<SizeType>(brep_json_format.size()));
                key("version");
                w_.Uint(brep_json_version);
                key("real");
                w_.String(real_name<T>);
                key("unit_scale");
                real(m_.unitScale());
                key("next_part");
                uint(m_.nextPart());
                key("geometry");
                w_.StartArray();
                for (std::size_t i = 0; i < nodes_.size(); ++i)
                {
                    ctx_ = std::format("geometry {}", i);
                    write_node(nodes_[i]);
                }
                w_.EndArray();
                table<brep::ShapeType::Vertex>();
                table<brep::ShapeType::Edge>();
                table<brep::ShapeType::Wire>();
                table<brep::ShapeType::Face>();
                table<brep::ShapeType::Shell>();
                table<brep::ShapeType::Solid>();
                table<brep::ShapeType::Compound>();
                w_.EndObject();
                if (!w_.IsComplete())
                    throw BrepWriteError{"the JSON text is not complete"};
            }
        };

        // =====================================================================
        // Reading
        // =====================================================================

        /// A place in the document, for the message of an error: a chain of members and indices up to the root.
        struct At
        {
            const At *up{};
            const char *key{}; ///< member name, or nullptr for an element of an array
            std::size_t index{};

            [[nodiscard]] At operator/(const char *k) const { return At{this, k, 0}; }
            template <std::integral I> // a template: 0 would otherwise convert to a const char * as well
            [[nodiscard]] At operator/(I i) const
            {
                return At{this, nullptr, static_cast<std::size_t>(i)};
            }

            /// The JSON pointer, as "/edges/12/u/0".
            [[nodiscard]] std::string pointer() const
            {
                std::vector<const At *> chain;
                for (auto p = this; p && p->up; p = p->up)
                    chain.push_back(p);
                std::string s;
                for (auto it = chain.rbegin(); it != chain.rend(); ++it)
                {
                    s += '/';
                    s += (*it)->key ? std::string{(*it)->key} : std::to_string((*it)->index);
                }
                return s;
            }
        };

        template <typename T>
        class Reader
        {
            using Geom = std::variant<std::shared_ptr<Curve<T, 2>>, std::shared_ptr<Curve<T, 3>>, std::shared_ptr<Surface<T, 2>>, std::shared_ptr<Surface<T, 3>>>;

            std::string_view source_;
            brep::Model<T> m_;
            std::vector<Geom> geom_;
            std::array<std::size_t, 7> sizes_{};
            std::uint64_t next_part_{1};

            [[noreturn]] void fail(const At &at, std::string_view why, BrepReadErrc code = BrepReadErrc::Invalid) const
            {
                throw BrepReadError{source_, code, at.pointer(), why};
            }

            // ---- values ----------------------------------------------------------------------

            struct Field
            {
                const char *key;
                bool required;
            };

            /// The members of an object, each one among `keys` and given once, the required ones
            /// present: their values in the order of `keys` (nullptr for one absent).
            template <std::size_t N>
            auto fields(const Value &v, const At &at, const std::array<Field, N> &keys) const -> std::array<const Value *, N>
            {
                if (!v.IsObject())
                    fail(at, "an object expected");
                std::array<const Value *, N> r{};
                for (auto it = v.MemberBegin(); it != v.MemberEnd(); ++it)
                {
                    const std::string_view name{it->name.GetString(), it->name.GetStringLength()};
                    std::size_t k = 0;
                    while (k < N && name != keys[k].key)
                        ++k;
                    if (k == N)
                        fail(at, std::format("unknown member \"{}\"", name));
                    if (r[k])
                        fail(at, std::format("member \"{}\" given twice", name));
                    r[k] = &it->value;
                }
                for (std::size_t k = 0; k < N; ++k)
                    if (keys[k].required && !r[k])
                        fail(at, std::format("member \"{}\" missing", keys[k].key));
                return r;
            }

            T real(const Value &v, const At &at) const
            {
                if (!v.IsNumber())
                    fail(at, "a number expected");
                const double d = v.GetDouble();
                if constexpr (std::same_as<T, float>)
                {
                    const auto f = static_cast<float>(d);
                    if (!same_bits(static_cast<double>(f), d) && !std::isnan(d))
                        fail(at, "a number of a float expected");
                    return f;
                }
                else
                    return d;
            }
            bool boolean(const Value &v, const At &at) const
            {
                if (!v.IsBool())
                    fail(at, "true or false expected");
                return v.GetBool();
            }
            std::uint64_t count(const Value &v, const At &at) const
            {
                if (!v.IsUint64())
                    fail(at, "an unsigned integer expected");
                return v.GetUint64();
            }
            const Value &array(const Value &v, const At &at) const
            {
                if (!v.IsArray())
                    fail(at, "an array expected");
                return v;
            }
            const Value &array(const Value &v, const At &at, SizeType n) const
            {
                if (!v.IsArray() || v.Size() != n)
                    fail(at, std::format("an array of {} expected", n));
                return v;
            }
            template <std::size_t N>
            auto reals(const Value &v, const At &at) const -> std::array<T, N>
            {
                array(v, at, N);
                std::array<T, N> r;
                for (std::size_t i = 0; i < N; ++i)
                    r[i] = real(v[static_cast<SizeType>(i)], at / i);
                return r;
            }
            auto reals(const Value &v, const At &at) const -> std::vector<T>
            {
                array(v, at);
                std::vector<T> r(v.Size());
                for (SizeType i = 0; i < v.Size(); ++i)
                    r[i] = real(v[i], at / i);
                return r;
            }
            template <std::size_t N>
            auto points(const Value &v, const At &at) const -> std::vector<std::array<T, N>>
            {
                array(v, at);
                std::vector<std::array<T, N>> r(v.Size());
                for (SizeType i = 0; i < v.Size(); ++i)
                    r[i] = reals<N>(v[i], at / i);
                return r;
            }
            Frame3<T> frame(const Value &v, const At &at) const
            {
                array(v, at, 4);
                Frame3<T> f;
                f.o = reals<3>(v[0], at / 0);
                f.x = reals<3>(v[1], at / 1);
                f.y = reals<3>(v[2], at / 2);
                f.z = reals<3>(v[3], at / 3);
                return f;
            }
            brep::Placement<T> placement(const Value &v, const At &at) const
            {
                const auto f = fields(v, at, std::array<Field, 2>{{{"r", true}, {"t", true}}});
                brep::Placement<T> p;
                const At ar = at / "r";
                array(*f[0], ar, 3);
                for (SizeType i = 0; i < 3; ++i)
                    p.r[i] = reals<3>((*f[0])[i], ar / i);
                p.t = reals<3>(*f[1], at / "t");
                return p;
            }
            brep::Orientation orient(const Value &v, const At &at) const
            {
                if (v.IsString() && v.GetStringLength() == 1)
                {
                    if (*v.GetString() == 'F')
                        return brep::Orientation::Forward;
                    if (*v.GetString() == 'R')
                        return brep::Orientation::Reversed;
                }
                fail(at, "\"F\" or \"R\" expected");
            }
            /// An id stored in an entity: null (invalid) or an index below the size of its table.
            template <brep::ShapeType K>
            brep::Id<K> id(const Value &v, const At &at) const
            {
                if (v.IsNull())
                    return {};
                const auto n = sizes_[std::to_underlying(K)];
                if (!v.IsUint() || v.GetUint() >= n)
                    fail(at, std::format("a {} index below {} or null expected", type_names[std::to_underlying(K)], n));
                return brep::Id<K>{v.GetUint()};
            }
            brep::ShapeId shape_id(const Value &v, const At &at) const
            {
                array(v, at, 2);
                const auto &t = v[0];
                if (t.IsString())
                {
                    const std::string_view name{t.GetString(), t.GetStringLength()};
                    const At ai = at / 1;
                    for (std::size_t k = 0; k < 7; ++k)
                        if (name == type_names[k])
                            switch (static_cast<brep::ShapeType>(k))
                            {
                            case brep::ShapeType::Vertex: return id<brep::ShapeType::Vertex>(v[1], ai);
                            case brep::ShapeType::Edge: return id<brep::ShapeType::Edge>(v[1], ai);
                            case brep::ShapeType::Wire: return id<brep::ShapeType::Wire>(v[1], ai);
                            case brep::ShapeType::Face: return id<brep::ShapeType::Face>(v[1], ai);
                            case brep::ShapeType::Shell: return id<brep::ShapeType::Shell>(v[1], ai);
                            case brep::ShapeType::Solid: return id<brep::ShapeType::Solid>(v[1], ai);
                            case brep::ShapeType::Compound: return id<brep::ShapeType::Compound>(v[1], ai);
                            }
                }
                fail(at / 0, "a type of shape expected (\"vertex\", \"edge\", \"wire\", \"face\", \"shell\", \"solid\" or \"compound\")");
            }

            // ---- geometry --------------------------------------------------------------------

            /// A geometry of the pool read so far, of class P (std::shared_ptr<Curve<T, 3>>...): null
            /// only if `nullable`.
            template <typename P>
            P geom(const Value &v, const At &at, bool nullable = false) const
            {
                if (nullable && v.IsNull())
                    return nullptr;
                if (!v.IsUint() || v.GetUint() >= geom_.size())
                    fail(at, std::format("an index of a geometry below {} expected{}", geom_.size(), nullable ? ", or null" : ""));
                const auto &g = geom_[v.GetUint()];
                if (auto p = std::get_if<P>(&g))
                    return *p;
                static constexpr Kind want = std::is_same_v<P, std::shared_ptr<Curve<T, 2>>>   ? Kind::Curve2
                                             : std::is_same_v<P, std::shared_ptr<Curve<T, 3>>> ? Kind::Curve3
                                             : std::is_same_v<P, std::shared_ptr<Surface<T, 2>>> ? Kind::Surface2
                                                                                                   : Kind::Surface3;
                fail(at, std::format("{} expected, geometry {} is {}", kind_name(want), v.GetUint(), kind_name(static_cast<Kind>(g.index()))));
            }
            template <std::size_t D>
            auto curve(const Value &v, const At &at) const { return geom<std::shared_ptr<Curve<T, D>>>(v, at); }
            template <std::size_t D>
            auto surface(const Value &v, const At &at) const { return geom<std::shared_ptr<Surface<T, D>>>(v, at); }

            /// The knots, degree and poles of a B-spline curve of N coordinates per pole, the bounds
            /// put back if given.
            template <typename C, std::size_t N>
            auto bspline_curve(const At &at, const Value &degree, const Value &knots, const Value &poles, const Value *bounds) const
                -> std::shared_ptr<C>
            {
                const auto p = count(degree, at / "degree");
                auto k = reals(knots, at / "knots");
                auto q = points<N>(poles, at / "poles");
                if (!valid_knots(k, q.size(), static_cast<std::size_t>(p)))
                    fail(at / "knots", std::format("{} knots in order, finite, for {} poles of degree {} expected", q.size() + p + 1, q.size(), p));
                auto c = std::make_shared<C>(q, k, static_cast<std::size_t>(p));
                if (bounds)
                {
                    const auto b = reals<2>(*bounds, at / "bounds");
                    c->trim(b[0], b[1], false);
                }
                return c;
            }
            BSCfunction<T> function(const Value &v, const At &at) const
            {
                const auto f = fields(v, at, std::array<Field, 4>{{{"degree", true}, {"knots", true}, {"poles", true}, {"bounds", false}}});
                return BSCfunction<T>{*bspline_curve<BSCurve<T, 1>, 1>(at, *f[0], *f[1], *f[2], f[3])};
            }

            template <std::size_t D>
            auto read_curve(std::string_view type, const Value &v, const At &at) -> std::shared_ptr<Curve<T, D>>
            {
                if (type == "bspline_curve" || type == "nurbs_curve")
                {
                    const auto f = fields(v, at, std::array<Field, 6>{{{"type", true}, {"dim", true}, {"degree", true}, {"knots", true}, {"poles", true}, {"bounds", false}}});
                    if (type == "bspline_curve")
                        return bspline_curve<BSCurve<T, D>, D>(at, *f[2], *f[3], *f[4], f[5]);
                    return bspline_curve<BSCurveRational<T, D>, D + 1>(at, *f[2], *f[3], *f[4], f[5]);
                }
                if (type == "trimmed_curve")
                {
                    const auto f = fields(v, at, std::array<Field, 4>{{{"type", true}, {"dim", true}, {"basis", true}, {"range", true}}});
                    auto b = curve<D>(*f[2], at / "basis");
                    const auto r = reals<2>(*f[3], at / "range");
                    return std::make_shared<CurveTrimmed<T, D>>(b, r[0], r[1]);
                }
                if (type == "composite_curve")
                {
                    const auto f = fields(v, at, std::array<Field, 3>{{{"type", true}, {"dim", true}, {"curves", true}}});
                    const At ac = at / "curves";
                    const auto &a = array(*f[2], ac);
                    if (a.Empty())
                        fail(ac, "at least one curve expected");
                    std::vector<std::shared_ptr<Curve<T, D>>> cs(a.Size());
                    for (SizeType i = 0; i < a.Size(); ++i)
                        cs[i] = curve<D>(a[i], ac / i);
                    return std::make_shared<CurveComposite<T, D>>(cs);
                }
                if (type == "reversed_curve")
                {
                    const auto f = fields(v, at, std::array<Field, 3>{{{"type", true}, {"dim", true}, {"basis", true}}});
                    return std::make_shared<CurveReversed<T, D>>(curve<D>(*f[2], at / "basis"));
                }
                if (type == "reparametrized_curve")
                {
                    const auto f = fields(v, at, std::array<Field, 4>{{{"type", true}, {"dim", true}, {"basis", true}, {"function", true}}});
                    return std::make_shared<CurveReparametrized<T, D>>(curve<D>(*f[2], at / "basis"), function(*f[3], at / "function"));
                }
                if (type == "extended_curve")
                {
                    const auto f = fields(v, at, std::array<Field, 4>{{{"type", true}, {"dim", true}, {"basis", true}, {"clamped", true}}});
                    return std::make_shared<CurveExtended<T, D>>(curve<D>(*f[2], at / "basis"), boolean(*f[3], at / "clamped"));
                }
                if (type == "curve_on_surface")
                {
                    const auto f = fields(v, at, std::array<Field, 4>{{{"type", true}, {"dim", true}, {"pcurve", true}, {"surface", true}}});
                    return std::make_shared<CurveOnSurface<T, D>>(curve<2>(*f[2], at / "pcurve"), surface<D>(*f[3], at / "surface"));
                }
                if (type == "offset_curve")
                {
                    if constexpr (D == 2)
                    {
                        const auto f = fields(v, at, std::array<Field, 4>{{{"type", true}, {"dim", true}, {"basis", true}, {"function", true}}});
                        return std::make_shared<CurveOffset2D<T, BSCfunction<T>>>(curve<2>(*f[2], at / "basis"), function(*f[3], at / "function"));
                    }
                    else
                    {
                        const auto f =
                            fields(v, at, std::array<Field, 5>{{{"type", true}, {"dim", true}, {"basis", true}, {"function", true}, {"direction", true}}});
                        const auto d = reals<3>(*f[4], at / "direction");
                        auto c = std::make_shared<CurveOffset3D<T, BSCfunction<T>>>(curve<3>(*f[2], at / "basis"), function(*f[3], at / "function"), d);
                        if (!same_bits(c->direction(), d))
                            fail(at / "direction", "a direction that its normalization gives back expected");
                        return c;
                    }
                }
                if constexpr (D == 3)
                {
                    if (type == "line")
                    {
                        const auto f = fields(v, at, std::array<Field, 5>{{{"type", true}, {"dim", true}, {"origin", true}, {"vector", true}, {"range", true}}});
                        const auto r = reals<2>(*f[4], at / "range");
                        return std::make_shared<LineCurve<T>>(reals<3>(*f[2], at / "origin"), reals<3>(*f[3], at / "vector"), r[0], r[1]);
                    }
                    if (type == "ellipse")
                    {
                        const auto f = fields(v, at, std::array<Field, 5>{{{"type", true}, {"dim", true}, {"radii", true}, {"frame", true}, {"range", true}}});
                        const auto r = reals<2>(*f[2], at / "radii");
                        const auto t = reals<2>(*f[4], at / "range");
                        return std::make_shared<EllipseCurve<T>>(r[0], r[1], frame(*f[3], at / "frame"), t[0], t[1]);
                    }
                    if (type == "moved_curve")
                    {
                        const auto f = fields(v, at, std::array<Field, 4>{{{"type", true}, {"dim", true}, {"basis", true}, {"placement", true}}});
                        auto b = curve<3>(*f[2], at / "basis");
                        if (typeid(*b) == typeid(brep::MovedCurve<T>)) // the constructor would compose the placements
                            fail(at / "basis", "a moved curve cannot be the basis of another one");
                        return std::make_shared<brep::MovedCurve<T>>(b, placement(*f[3], at / "placement"));
                    }
                }
                fail(at / "type", std::format("no {}D curve of type \"{}\"", D, type));
            }

            auto elementary(const Value &frame_v, const Value &box_v, const At &at, auto &&make) const -> std::shared_ptr<Surface<T, 3>>
            {
                const auto f = frame(frame_v, at / "frame");
                const auto box = reals<4>(box_v, at / "box");
                // made on a box every primitive accepts, then given its own: `with_box`, as the
                // operations that move a seam, checks only that the box is finite and not empty
                // (a sphere moved there may run beyond the poles by an ulp)
                std::shared_ptr<ElementarySurface<T>> s = make(f, std::array<T, 4>{T(0), T(1), T(0), T(1)});
                return s->with_box(box);
            }

            template <std::size_t D>
            auto read_surface(std::string_view type, const Value &v, const At &at) -> std::shared_ptr<Surface<T, D>>
            {
                if (type == "bspline_surface" || type == "nurbs_surface")
                {
                    const auto f = fields(v, at, std::array<Field, 6>{{{"type", true}, {"dim", true}, {"degree", true}, {"knots", true}, {"poles", true}, {"bounds", false}}});
                    const At ad = at / "degree", ak = at / "knots";
                    array(*f[2], ad, 2);
                    array(*f[3], ak, 2);
                    const auto p = count((*f[2])[0], ad / 0), q = count((*f[2])[1], ad / 1);
                    auto ku = reals((*f[3])[0], ak / 0), kv = reals((*f[3])[1], ak / 1);
                    const auto nu = ku.size() > p + 1 ? ku.size() - p - 1 : 0, nv = kv.size() > q + 1 ? kv.size() - q - 1 : 0;
                    if (!valid_knots(ku, nu, static_cast<std::size_t>(p)))
                        fail(ak / 0, std::format("knots in order, finite, more than twice the degree {} plus one, expected", p));
                    if (!valid_knots(kv, nv, static_cast<std::size_t>(q)))
                        fail(ak / 1, std::format("knots in order, finite, more than twice the degree {} plus one, expected", q));
                    auto make = [&]<bool R>() -> std::shared_ptr<Surface<T, D>> {
                        using S = std::conditional_t<R, BSSurfaceRational<T, D>, BSSurface<T, D>>;
                        auto poles = points<D + R>(*f[4], at / "poles");
                        if (poles.size() != nu * nv)
                            fail(at / "poles", std::format("{} x {} poles expected", nu, nv));
                        auto s = std::make_shared<S>(poles, ku, kv, static_cast<std::size_t>(p), static_cast<std::size_t>(q));
                        if (f[5])
                        {
                            const auto b = reals<4>(*f[5], at / "bounds");
                            if (!same_bits(b[2], kv.front()) || !same_bits(b[3], kv.back()))
                                fail(at / "bounds", "V bounds equal to the end knots expected");
                            if (!trim_u_noop(ku, b[0], b[1]))
                                fail(at / "bounds", "U bounds within knot_eps of the end knots expected");
                            s->trimU(b[0], b[1]);
                        }
                        return s;
                    };
                    return type == "bspline_surface" ? make.template operator()<false>() : make.template operator()<true>();
                }
                if constexpr (D == 3)
                {
                    if (type == "plane")
                    {
                        const auto f = fields(v, at, std::array<Field, 4>{{{"type", true}, {"dim", true}, {"frame", true}, {"box", true}}});
                        return elementary(*f[2], *f[3], at, [](const Frame3<T> &fr, std::array<T, 4> b) { return std::make_shared<Plane<T>>(fr, b); });
                    }
                    if (type == "cylinder")
                    {
                        const auto f = fields(v, at, std::array<Field, 5>{{{"type", true}, {"dim", true}, {"radius", true}, {"frame", true}, {"box", true}}});
                        const T r = real(*f[2], at / "radius");
                        return elementary(*f[3], *f[4], at, [r](const Frame3<T> &fr, std::array<T, 4> b) {
                            return std::make_shared<CylindricalSurface<T>>(r, fr, b);
                        });
                    }
                    if (type == "cone")
                    {
                        const auto f = fields(
                            v, at, std::array<Field, 6>{{{"type", true}, {"dim", true}, {"radius", true}, {"semi_angle", true}, {"frame", true}, {"box", true}}});
                        const T r = real(*f[2], at / "radius"), a = real(*f[3], at / "semi_angle");
                        return elementary(*f[4], *f[5], at, [r, a](const Frame3<T> &fr, std::array<T, 4> b) {
                            return std::make_shared<ConicalSurface<T>>(r, a, fr, b);
                        });
                    }
                    if (type == "sphere")
                    {
                        const auto f = fields(v, at, std::array<Field, 5>{{{"type", true}, {"dim", true}, {"radius", true}, {"frame", true}, {"box", true}}});
                        const T r = real(*f[2], at / "radius");
                        return elementary(*f[3], *f[4], at, [r](const Frame3<T> &fr, std::array<T, 4> b) {
                            return std::make_shared<SphericalSurface<T>>(r, fr, b);
                        });
                    }
                    if (type == "torus")
                    {
                        const auto f = fields(
                            v, at,
                            std::array<Field, 6>{{{"type", true}, {"dim", true}, {"major_radius", true}, {"minor_radius", true}, {"frame", true}, {"box", true}}});
                        const T a = real(*f[2], at / "major_radius"), b = real(*f[3], at / "minor_radius");
                        return elementary(*f[4], *f[5], at, [a, b](const Frame3<T> &fr, std::array<T, 4> box) {
                            return std::make_shared<ToroidalSurface<T>>(a, b, fr, box);
                        });
                    }
                    if (type == "revolution_surface")
                    {
                        const auto f = fields(v, at, std::array<Field, 5>{{{"type", true}, {"dim", true}, {"curve", true}, {"axis", true}, {"angles", true}}});
                        auto c = curve<2>(*f[2], at / "curve");
                        const At aa = at / "axis";
                        array(*f[3], aa, 3);
                        ax2<T, 3> ax{reals<3>((*f[3])[0], aa / 0), reals<3>((*f[3])[1], aa / 1), reals<3>((*f[3])[2], aa / 2)};
                        const auto a = reals<2>(*f[4], at / "angles");
                        return std::make_shared<SurfaceOfRevolution<T>>(c, ax, a[0], a[1]);
                    }
                    if (type == "offset_surface")
                    {
                        const auto f = fields(v, at, std::array<Field, 4>{{{"type", true}, {"dim", true}, {"basis", true}, {"distance", true}}});
                        return std::make_shared<OffsetSurface<T>>(surface<3>(*f[2], at / "basis"), real(*f[3], at / "distance"));
                    }
                    if (type == "moved_surface")
                    {
                        const auto f = fields(v, at, std::array<Field, 4>{{{"type", true}, {"dim", true}, {"basis", true}, {"placement", true}}});
                        auto b = surface<3>(*f[2], at / "basis");
                        if (typeid(*b) == typeid(brep::MovedSurface<T>)) // the constructor would compose the placements
                            fail(at / "basis", "a moved surface cannot be the basis of another one");
                        return std::make_shared<brep::MovedSurface<T>>(b, placement(*f[3], at / "placement"));
                    }
                    if (type == "seamless_surface")
                    {
                        const auto f = fields(v, at, std::array<Field, 4>{{{"type", true}, {"dim", true}, {"basis", true}, {"periods", true}}});
                        return std::make_shared<brep::SeamlessSurface<T>>(surface<3>(*f[2], at / "basis"), reals<2>(*f[3], at / "periods"));
                    }
                }
                fail(at / "type", std::format("no {}D surface of type \"{}\"", D, type));
            }

            static bool is_surface_type(std::string_view t)
            {
                for (const char *s : {"bspline_surface", "nurbs_surface", "plane", "cylinder", "cone", "sphere", "torus", "revolution_surface", "offset_surface",
                                      "moved_surface", "seamless_surface"})
                    if (t == s)
                        return true;
                return false;
            }

            Geom read_geom(const Value &v, const At &at)
            {
                if (!v.IsObject())
                    fail(at, "an object expected");
                const auto t = v.FindMember("type");
                if (t == v.MemberEnd() || !t->value.IsString())
                    fail(at / "type", "the type of the geometry, a string, expected");
                const auto d = v.FindMember("dim");
                if (d == v.MemberEnd() || !d->value.IsUint() || (d->value.GetUint() != 2 && d->value.GetUint() != 3))
                    fail(at / "dim", "2 or 3 expected");
                const std::string_view type{t->value.GetString(), t->value.GetStringLength()};
                try
                {
                    if (is_surface_type(type))
                    {
                        if (d->value.GetUint() == 2)
                            return read_surface<2>(type, v, at);
                        return read_surface<3>(type, v, at);
                    }
                    if (d->value.GetUint() == 2)
                        return read_curve<2>(type, v, at);
                    return read_curve<3>(type, v, at);
                }
                catch (const std::invalid_argument &e) // a constructor of gbs refused the values
                {
                    fail(at, e.what());
                }
                catch (const std::length_error &e)
                {
                    fail(at, e.what());
                }
            }

            // ---- topology --------------------------------------------------------------------

            brep::Vertex<T> vertex(const Value &v, const At &at, const Value *(&attrs)[3]) const
            {
                const auto f = fields(v, at, std::array<Field, 5>{{{"p", true}, {"tol", true}, {"name", false}, {"external_id", false}, {"instance", false}}});
                attrs[0] = f[2], attrs[1] = f[3], attrs[2] = f[4];
                return brep::Vertex<T>{reals<3>(*f[0], at / "p"), real(*f[1], at / "tol")};
            }
            brep::Edge<T> edge(const Value &v, const At &at, const Value *(&attrs)[3]) const
            {
                const auto f = fields(v, at, std::array<Field, 9>{{{"curve", true}, {"u", true}, {"v", true}, {"tol", true}, {"degenerate", true},
                                                                   {"same_parameter", true}, {"name", false}, {"external_id", false}, {"instance", false}}});
                attrs[0] = f[6], attrs[1] = f[7], attrs[2] = f[8];
                brep::Edge<T> e;
                e.curve = geom<std::shared_ptr<Curve<T, 3>>>(*f[0], at / "curve", true);
                const auto u = reals<2>(*f[1], at / "u");
                e.u1 = u[0];
                e.u2 = u[1];
                const At av = at / "v";
                array(*f[2], av, 2);
                e.v1 = id<brep::ShapeType::Vertex>((*f[2])[0], av / 0);
                e.v2 = id<brep::ShapeType::Vertex>((*f[2])[1], av / 1);
                e.tol = real(*f[3], at / "tol");
                e.degenerate = boolean(*f[4], at / "degenerate");
                e.same_parameter = boolean(*f[5], at / "same_parameter");
                return e;
            }
            brep::Wire<T> wire(const Value &v, const At &at, const Value *(&attrs)[3]) const
            {
                const auto f = fields(v, at, std::array<Field, 5>{{{"coedges", true}, {"closed", true}, {"name", false}, {"external_id", false}, {"instance", false}}});
                attrs[0] = f[2], attrs[1] = f[3], attrs[2] = f[4];
                brep::Wire<T> w;
                const At ac = at / "coedges";
                const auto &a = array(*f[0], ac);
                w.coedges.resize(a.Size());
                for (SizeType i = 0; i < a.Size(); ++i)
                {
                    const At ai = ac / i;
                    array(a[i], ai, 3);
                    auto &ce = w.coedges[i];
                    ce.edge = id<brep::ShapeType::Edge>(a[i][0], ai / 0);
                    ce.orient = orient(a[i][1], ai / 1);
                    ce.pcurve = geom<std::shared_ptr<Curve<T, 2>>>(a[i][2], ai / 2, true);
                }
                w.closed = boolean(*f[1], at / "closed");
                return w;
            }
            brep::Face<T> face(const Value &v, const At &at, const Value *(&attrs)[3]) const
            {
                const auto f = fields(v, at, std::array<Field, 7>{{{"surface", true}, {"wires", true}, {"tol", true}, {"natural_bounds", true},
                                                                   {"name", false}, {"external_id", false}, {"instance", false}}});
                attrs[0] = f[4], attrs[1] = f[5], attrs[2] = f[6];
                brep::Face<T> fc;
                fc.surface = geom<std::shared_ptr<Surface<T, 3>>>(*f[0], at / "surface", true);
                const At aw = at / "wires";
                const auto &a = array(*f[1], aw);
                fc.wires.resize(a.Size());
                for (SizeType i = 0; i < a.Size(); ++i)
                    fc.wires[i] = id<brep::ShapeType::Wire>(a[i], aw / i);
                fc.tol = real(*f[2], at / "tol");
                fc.natural_bounds = boolean(*f[3], at / "natural_bounds");
                return fc;
            }
            brep::Shell shell(const Value &v, const At &at, const Value *(&attrs)[3]) const
            {
                const auto f = fields(v, at, std::array<Field, 5>{{{"faces", true}, {"closed", true}, {"name", false}, {"external_id", false}, {"instance", false}}});
                attrs[0] = f[2], attrs[1] = f[3], attrs[2] = f[4];
                brep::Shell s;
                const At af = at / "faces";
                const auto &a = array(*f[0], af);
                s.faces.resize(a.Size());
                for (SizeType i = 0; i < a.Size(); ++i)
                {
                    const At ai = af / i;
                    array(a[i], ai, 2);
                    s.faces[i].face = id<brep::ShapeType::Face>(a[i][0], ai / 0);
                    s.faces[i].orient = orient(a[i][1], ai / 1);
                }
                s.closed = boolean(*f[1], at / "closed");
                return s;
            }
            brep::Solid solid(const Value &v, const At &at, const Value *(&attrs)[3]) const
            {
                const auto f = fields(v, at, std::array<Field, 5>{{{"outer", true}, {"voids", true}, {"name", false}, {"external_id", false}, {"instance", false}}});
                attrs[0] = f[2], attrs[1] = f[3], attrs[2] = f[4];
                brep::Solid s;
                s.outer = id<brep::ShapeType::Shell>(*f[0], at / "outer");
                const At av = at / "voids";
                const auto &a = array(*f[1], av);
                s.voids.resize(a.Size());
                for (SizeType i = 0; i < a.Size(); ++i)
                    s.voids[i] = id<brep::ShapeType::Shell>(a[i], av / i);
                return s;
            }
            brep::Compound compound(const Value &v, const At &at, const Value *(&attrs)[3]) const
            {
                const auto f = fields(v, at, std::array<Field, 4>{{{"shapes", true}, {"name", false}, {"external_id", false}, {"instance", false}}});
                attrs[0] = f[1], attrs[1] = f[2], attrs[2] = f[3];
                brep::Compound c;
                const At as_ = at / "shapes";
                const auto &a = array(*f[0], as_);
                c.shapes.reserve(a.Size());
                for (SizeType i = 0; i < a.Size(); ++i)
                    c.shapes.push_back(shape_id(a[i], as_ / i));
                return c;
            }

            /// name, external id, instance of the entity just added
            void attributes(const brep::ShapeId &s, const Value *const (&attrs)[3], const At &at)
            {
                if (attrs[0])
                {
                    const At an = at / "name";
                    if (!attrs[0]->IsString() || attrs[0]->GetStringLength() == 0)
                        fail(an, "a name, a string not empty, expected");
                    m_.setName(s, std::string{attrs[0]->GetString(), attrs[0]->GetStringLength()});
                }
                if (attrs[1])
                {
                    if (!attrs[1]->IsInt64())
                        fail(at / "external_id", "an integer expected");
                    m_.setExternalId(s, attrs[1]->GetInt64());
                }
                if (attrs[2])
                {
                    const At ai = at / "instance";
                    const auto f = fields(*attrs[2], ai, std::array<Field, 2>{{{"part", true}, {"placement", true}}});
                    const auto part = count(*f[0], ai / "part");
                    if (part == 0 || part >= next_part_)
                        fail(ai / "part", std::format("a part key from 1 to {} expected (below next_part)", next_part_ - 1));
                    m_.setInstance(s, brep::InstanceOf<T>{part, placement(*f[1], ai / "placement")});
                }
            }

            template <brep::ShapeType K>
            void table(const Value &a, const At &at)
            {
                for (SizeType i = 0; i < a.Size(); ++i)
                {
                    const At ai = at / i;
                    const Value *attrs[3]{};
                    const auto &v = a[i];
                    brep::Id<K> h;
                    if constexpr (K == brep::ShapeType::Vertex) h = v.IsNull() ? m_.add(brep::Vertex<T>{}) : m_.add(vertex(v, ai, attrs));
                    else if constexpr (K == brep::ShapeType::Edge) h = v.IsNull() ? m_.add(brep::Edge<T>{}) : m_.add(edge(v, ai, attrs));
                    else if constexpr (K == brep::ShapeType::Wire) h = v.IsNull() ? m_.add(brep::Wire<T>{}) : m_.add(wire(v, ai, attrs));
                    else if constexpr (K == brep::ShapeType::Face) h = v.IsNull() ? m_.add(brep::Face<T>{}) : m_.add(face(v, ai, attrs));
                    else if constexpr (K == brep::ShapeType::Shell) h = v.IsNull() ? m_.add(brep::Shell{}) : m_.add(shell(v, ai, attrs));
                    else if constexpr (K == brep::ShapeType::Solid) h = v.IsNull() ? m_.add(brep::Solid{}) : m_.add(solid(v, ai, attrs));
                    else h = v.IsNull() ? m_.add(brep::Compound{}) : m_.add(compound(v, ai, attrs));
                    if (v.IsNull())
                        m_.erase(h); // dead, as it was
                    else
                        attributes(h, attrs, ai);
                }
            }

        public:
            explicit Reader(std::string_view source) : source_{source} {}

            brep::Model<T> run(const Value &root)
            {
                const At at{};
                if (!root.IsObject())
                    fail(at, "not a BREP model of gbs: a JSON object expected", BrepReadErrc::NotBrep);
                // format and version first: a newer version may hold what this code does not know
                const auto fmt = root.FindMember("format");
                if (fmt == root.MemberEnd() || !fmt->value.IsString() ||
                    std::string_view{fmt->value.GetString(), fmt->value.GetStringLength()} != brep_json_format)
                    fail(at / "format", std::format("not a BREP model of gbs: \"format\": \"{}\" expected", brep_json_format), BrepReadErrc::NotBrep);
                const auto ver = root.FindMember("version");
                if (ver == root.MemberEnd() || !ver->value.IsUint() || ver->value.GetUint() == 0)
                    fail(at / "version", "the version of the format, an integer from 1, expected");
                if (ver->value.GetUint() > brep_json_version)
                    fail(at / "version",
                         std::format("version {} of the format, newer than this reader (version {}): a newer gbs is needed", ver->value.GetUint(),
                                     brep_json_version),
                         BrepReadErrc::NewerVersion);
                const auto f = fields(root, at,
                                      std::array<Field, 13>{{{"format", true},
                                                             {"version", true},
                                                             {"real", true},
                                                             {"unit_scale", true},
                                                             {"next_part", true},
                                                             {"geometry", true},
                                                             {"vertices", true},
                                                             {"edges", true},
                                                             {"wires", true},
                                                             {"faces", true},
                                                             {"shells", true},
                                                             {"solids", true},
                                                             {"compounds", true}}});
                if (!f[2]->IsString() || std::string_view{f[2]->GetString(), f[2]->GetStringLength()} != real_name<T>)
                    fail(at / "real", std::format("\"{}\" expected: the model is read in {}", real_name<T>, real_name<T>));
                const T unit = real(*f[3], at / "unit_scale");
                if (!std::isfinite(unit) || !(unit > T(0)))
                    fail(at / "unit_scale", "a length of a unit finite and > 0 expected");
                m_.setUnitScale(unit);
                next_part_ = count(*f[4], at / "next_part");
                if (next_part_ == 0)
                    fail(at / "next_part", "a part key from 1 expected");
                m_.reservePartKeys(next_part_);
                for (std::size_t k = 0; k < 7; ++k)
                {
                    const At ak = at / table_keys[k];
                    const auto &a = array(*f[6 + k], ak);
                    if (a.Size() >= brep::npos)
                        fail(ak, "too many entities");
                    sizes_[k] = a.Size();
                }
                const At ag = at / "geometry";
                const auto &g = array(*f[5], ag);
                geom_.reserve(g.Size());
                for (SizeType i = 0; i < g.Size(); ++i)
                    geom_.push_back(read_geom(g[i], ag / i));
                table<brep::ShapeType::Vertex>(*f[6], at / table_keys[0]);
                table<brep::ShapeType::Edge>(*f[7], at / table_keys[1]);
                table<brep::ShapeType::Wire>(*f[8], at / table_keys[2]);
                table<brep::ShapeType::Face>(*f[9], at / table_keys[3]);
                table<brep::ShapeType::Shell>(*f[10], at / table_keys[4]);
                table<brep::ShapeType::Solid>(*f[11], at / table_keys[5]);
                table<brep::ShapeType::Compound>(*f[12], at / table_keys[6]);
                return std::move(m_);
            }
        };

        /// Parses the text, without recursion, numbers in full precision, NaN and infinities allowed.
        inline void parse(std::string_view text, rapidjson::Document &d, std::string_view source)
        {
            constexpr unsigned flags = rapidjson::kParseFullPrecisionFlag | rapidjson::kParseIterativeFlag | rapidjson::kParseNanAndInfFlag;
            d.Parse<flags>(text.data(), text.size());
            if (d.HasParseError())
            {
                const auto offset = std::min<std::size_t>(d.GetErrorOffset(), text.size());
                const auto before = text.substr(0, offset);
                const auto line_start = before.rfind('\n');
                const auto line = static_cast<std::size_t>(std::ranges::count(before, '\n')) + 1;
                const auto column = offset - (line_start == std::string_view::npos ? 0 : line_start + 1) + 1;
                throw BrepReadError{source, BrepReadErrc::Syntax, "",
                                    std::format("line {}, column {}: {}", line, column, rapidjson::GetParseError_En(d.GetParseError()))};
            }
            if (const auto nul = text.find('\0'); nul != std::string_view::npos)
                throw BrepReadError{source, BrepReadErrc::Syntax, "", std::format("NUL byte at offset {}", nul)};
        }

        template <typename T, typename W>
        void write(const brep::Model<T> &m, W &w)
        {
            Writer<T, W>{w, m}.run();
        }
    }

    // =========================================================================
    // API
    // =========================================================================

    /**
     * @brief The JSON text of a model, read back exactly by `read_brep`.
     * @throws BrepWriteError for a geometry of a class the format does not hold (an application's
     * curve, the old `Line` and `Circle`), an id stored in an entity beyond its table
     */
    template <brep_json_real T>
    [[nodiscard]] auto write_brep(const brep::Model<T> &m, const BrepWriteOptions &options = {}) -> std::string
    {
        rapidjson::StringBuffer buffer;
        using Enc = rapidjson::UTF8<>;
        if (options.pretty)
        {
            rapidjson::PrettyWriter<rapidjson::StringBuffer, Enc, Enc, rapidjson::CrtAllocator, rapidjson::kWriteNanAndInfFlag> w(buffer);
            w.SetIndent(' ', 1);
            w.SetFormatOptions(rapidjson::kFormatSingleLineArray);
            brep_json_detail::write(m, w);
        }
        else
        {
            rapidjson::Writer<rapidjson::StringBuffer, Enc, Enc, rapidjson::CrtAllocator, rapidjson::kWriteNanAndInfFlag> w(buffer);
            brep_json_detail::write(m, w);
        }
        return std::string{buffer.GetString(), buffer.GetSize()};
    }

    /**
     * @brief Writes the JSON text of a model to a stream.
     * @throws BrepWriteError as `write_brep`; std::ios_base::failure if the stream fails
     */
    template <brep_json_real T>
    void write_brep(const brep::Model<T> &m, std::ostream &os, const BrepWriteOptions &options = {})
    {
        const auto text = write_brep(m, options);
        os.write(text.data(), static_cast<std::streamsize>(text.size()));
        if (!os)
            throw std::ios_base::failure{"write_brep: cannot write the stream"};
    }

    /**
     * @brief Writes the JSON text of a model to a file, whole or not at all (`write_text_atomic`).
     * @throws BrepWriteError as `write_brep` (the file untouched); std::filesystem::filesystem_error
     * if the file cannot be written
     */
    template <brep_json_real T>
    void write_brep_file(const brep::Model<T> &m, const std::filesystem::path &file, const BrepWriteOptions &options = {})
    {
        write_text_atomic(file, write_brep(m, options));
    }

    /**
     * @brief The model of a JSON text written by `write_brep`: the same tables, ids, sharing and
     * bits.
     * @param source what the text is, for the message of an error (a file name)
     * @throws BrepReadError (code, where, reason) for a text that is not such a model
     */
    template <brep_json_real T = double>
    [[nodiscard]] auto read_brep(std::string_view text, std::string_view source = "BREP JSON") -> brep::Model<T>
    {
        rapidjson::Document d;
        brep_json_detail::parse(text, d, source);
        return brep_json_detail::Reader<T>{source}.run(d);
    }

    /**
     * @brief The model of a stream holding a JSON text written by `write_brep`, read to its end.
     * @throws BrepReadError as `read_brep`; std::ios_base::failure if the stream cannot be read
     */
    template <brep_json_real T = double>
    [[nodiscard]] auto read_brep(std::istream &is, std::string_view source = "BREP JSON stream") -> brep::Model<T>
    {
        const std::string text{std::istreambuf_iterator<char>{is}, std::istreambuf_iterator<char>{}};
        if (is.bad())
            throw std::ios_base::failure{"read_brep: cannot read the stream"};
        return read_brep<T>(text, source);
    }

    /**
     * @brief The model of a file written by `write_brep_file`.
     * @throws BrepReadError as `read_brep` (the source is the file name);
     * std::filesystem::filesystem_error if the file cannot be read
     */
    template <brep_json_real T = double>
    [[nodiscard]] auto read_brep_file(const std::filesystem::path &file) -> brep::Model<T>
    {
        return read_brep<T>(read_text_or_throw(file, "cannot read the BREP file"), utf8_string(file));
    }
}
