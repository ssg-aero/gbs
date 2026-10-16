#pragma once
/**
 * @file curve_pieces.h
 * @brief Pieces of a curve for subdivision algorithms: a box enclosing each piece, its
 * flatness, and the split of a piece in two halves.
 *
 * A B-spline (rational or not) is cut into Bézier segments: the box of the Cartesian poles
 * of a segment encloses it (convex hull property, positive weights), and a segment is halved
 * exactly by de Casteljau. Any other curve is sampled: its box is the box of the samples,
 * grown by a margin estimated from the deviation of the curve to the chords of the samples.
 * That box is not guaranteed; it is wide enough for curves whose curvature does not vary
 * abruptly between samples.
 */
#include <gbs/intersect/box.h>
#include <gbs/bscurve.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#ifdef GBS_USE_MODULES
    import vecop;
#endif

namespace gbs
{
    /// A piece of a curve on the parameter range [a, b].
    template <std::floating_point T, std::size_t dim>
    struct CurvePiece
    {
        T a{}, b{};
        /// Bézier piece: homogeneous poles (w x, w); sampled piece: (x, 1) at uniform parameters
        std::vector<std::array<T, dim + 1>> poles;
        bool bezier{};
        Box<T, dim> box;  ///< encloses the piece
        T flatness{};     ///< bound of the distance of the piece to the chord of its end points
        point<T, dim> first{}, last{};

        [[nodiscard]] T chord() const { return norm(last - first); }
        [[nodiscard]] T size() const { return box.diagonal(); }
    };

    namespace intersect_detail
    {
        template <std::floating_point T, std::size_t dim>
        auto cartesian(const std::array<T, dim + 1> &h) -> point<T, dim>
        {
            point<T, dim> p;
            for (std::size_t i{}; i < dim; ++i)
                p[i] = h[i] / h[dim];
            return p;
        }

        /// Distance of `p` to the segment [a, b] (to `a` when the segment is a point).
        template <std::floating_point T, std::size_t dim>
        auto distance_to_segment(const point<T, dim> &p, const point<T, dim> &a, const point<T, dim> &b) -> T
        {
            const auto ab = b - a;
            const auto l2 = sq_norm(ab);
            if (l2 <= T(0))
                return norm(p - a);
            const auto s = std::clamp((p - a) * ab / l2, T(0), T(1));
            return norm(p - (a + s * ab));
        }
    }

    /// Cuts a curve into pieces and halves them on demand.
    template <std::floating_point T, std::size_t dim>
    class CurvePieces
    {
        const Curve<T, dim> &m_crv;
        std::size_t m_samples;

        /**
         * Bézier segments of a clamped B-spline (The NURBS Book, algorithm A5.6): every interior
         * knot raised to the degree at once, in O(n p²), instead of one insertion at a time.
         */
        template <typename Pole>
        static auto decompose(const std::vector<T> &U, const std::vector<Pole> &Pw, std::size_t p) -> std::vector<std::vector<Pole>>
        {
            const std::size_t n = Pw.size() - 1, m = n + p + 1;
            std::vector<std::vector<Pole>> Q(1, std::vector<Pole>(p + 1));
            for (std::size_t i = 0; i <= p; ++i)
                Q[0][i] = Pw[i];
            std::vector<T> alphas(p + 1);
            std::size_t a = p, b = p + 1;
            while (b < m)
            {
                const std::size_t i0 = b;
                while (b < m && U[b + 1] == U[b])
                    ++b;
                const std::size_t mult = b - i0 + 1;
                auto &cur = Q.back();
                std::vector<Pole> next(p + 1);
                if (mult < p)
                {
                    const T numer = U[b] - U[a];
                    for (std::size_t j = p; j > mult; --j)
                        alphas[j - mult - 1] = numer / (U[a + j] - U[a]);
                    const std::size_t r = p - mult;
                    for (std::size_t j = 1; j <= r; ++j)
                    {
                        const std::size_t save = r - j, s = mult + j;
                        for (std::size_t k = p; k >= s; --k)
                        {
                            const T alpha = alphas[k - s];
                            for (std::size_t c = 0; c < cur[k].size(); ++c)
                                cur[k][c] = alpha * cur[k][c] + (T(1) - alpha) * cur[k - 1][c];
                        }
                        if (b < m)
                            next[save] = cur[p];
                    }
                }
                if (b < m)
                {
                    for (std::size_t i = p - mult; i <= p; ++i)
                        next[i] = Pw[b - p + i];
                    Q.push_back(std::move(next));
                    a = b;
                    ++b;
                }
            }
            return Q;
        }

        template <typename Bs>
        auto bezier_pieces(Bs c, T a, T b) -> std::vector<CurvePiece<T, dim>>
        {
            c.trim(a, b);
            const auto p = c.degree();
            const auto knots = c.knots();
            const auto &U = c.knotsFlats();
            const auto n_seg = knots.size() - 1;
            std::vector<CurvePiece<T, dim>> pieces;
            if (p < 1 || U.size() != c.poles().size() + p + 1 || U[p] != U.front() || U[U.size() - p - 1] != U.back())
                return pieces; // not clamped: the caller samples the curve
            const auto segments = decompose(U, c.poles(), p);
            if (segments.size() != n_seg)
                return pieces;
            pieces.reserve(n_seg);
            for (std::size_t s{}; s < n_seg; ++s)
            {
                CurvePiece<T, dim> pc{.a = knots[s], .b = knots[s + 1], .bezier = true};
                pc.poles.resize(p + 1);
                for (std::size_t k{}; k <= p; ++k)
                {
                    const auto &q = segments[s][k];
                    if constexpr (std::tuple_size_v<std::remove_cvref_t<decltype(q)>> == dim + 1)
                        pc.poles[k] = q;
                    else
                    {
                        for (std::size_t i{}; i < dim; ++i)
                            pc.poles[k][i] = q[i];
                        pc.poles[k][dim] = T(1);
                    }
                }
                finish_bezier(pc);
                pieces.push_back(std::move(pc));
            }
            return pieces;
        }

        static void finish_bezier(CurvePiece<T, dim> &pc)
        {
            pc.box = Box<T, dim>{};
            pc.first = intersect_detail::cartesian<T, dim>(pc.poles.front());
            pc.last = intersect_detail::cartesian<T, dim>(pc.poles.back());
            pc.flatness = T(0);
            for (const auto &h : pc.poles)
            {
                const auto x = intersect_detail::cartesian<T, dim>(h);
                pc.box.add(x);
                pc.flatness = std::max(pc.flatness, intersect_detail::distance_to_segment(x, pc.first, pc.last));
            }
        }

        auto sampled_piece(T a, T b) const -> CurvePiece<T, dim>
        {
            CurvePiece<T, dim> pc{.a = a, .b = b, .bezier = false};
            const auto n = m_samples;
            pc.poles.resize(n + 1);
            std::vector<point<T, dim>> x(n + 1);
            for (std::size_t k{}; k <= n; ++k)
            {
                x[k] = m_crv.value(k == n ? b : a + (b - a) * T(k) / T(n));
                for (std::size_t i{}; i < dim; ++i)
                    pc.poles[k][i] = x[k][i];
                pc.poles[k][dim] = T(1);
                pc.box.add(x[k]);
            }
            pc.first = x.front();
            pc.last = x.back();
            // margin: twice the largest deviation of the curve to the chords of the samples
            T margin{};
            for (std::size_t k{}; k < n; ++k)
            {
                const auto mid = m_crv.value(a + (b - a) * (T(k) + T(0.5)) / T(n));
                margin = std::max(margin, intersect_detail::distance_to_segment(mid, x[k], x[k + 1]));
            }
            margin *= T(2);
            pc.box.inflate(margin);
            for (const auto &p : x)
                pc.flatness = std::max(pc.flatness, intersect_detail::distance_to_segment(p, pc.first, pc.last));
            pc.flatness += margin;
            return pc;
        }

    public:
        /// @param samples number of intervals sampled on a piece of a curve which is not a B-spline
        explicit CurvePieces(const Curve<T, dim> &crv, std::size_t samples = 8) : m_crv{crv}, m_samples{std::max<std::size_t>(samples, 2)} {}

        [[nodiscard]] auto curve() const -> const Curve<T, dim> & { return m_crv; }

        /// The initial pieces of the curve on [a, b]: its Bézier segments, or one sampled piece.
        [[nodiscard]] auto pieces(T a, T b) -> std::vector<CurvePiece<T, dim>>
        {
            std::vector<CurvePiece<T, dim>> r;
            if (const auto *bs = dynamic_cast<const BSCurve<T, dim> *>(&m_crv))
                r = bezier_pieces(*bs, a, b);
            else if (const auto *bsr = dynamic_cast<const BSCurveRational<T, dim> *>(&m_crv))
                r = bezier_pieces(*bsr, a, b);
            if (r.empty())
                r.push_back(sampled_piece(a, b));
            return r;
        }

        /// The two halves of a piece, at the middle of its parameter range.
        [[nodiscard]] auto split(const CurvePiece<T, dim> &pc) const -> std::pair<CurvePiece<T, dim>, CurvePiece<T, dim>>
        {
            const auto m = (pc.a + pc.b) / T(2);
            if (!pc.bezier)
                return {sampled_piece(pc.a, m), sampled_piece(m, pc.b)};
            // de Casteljau at t = 1/2, in homogeneous coordinates
            const auto n = pc.poles.size();
            CurvePiece<T, dim> l{.a = pc.a, .b = m, .bezier = true}, r{.a = m, .b = pc.b, .bezier = true};
            l.poles.resize(n);
            r.poles.resize(n);
            auto w = pc.poles;
            for (std::size_t k{}; k < n; ++k)
            {
                l.poles[k] = w.front();
                r.poles[n - 1 - k] = w[n - 1 - k];
                for (std::size_t j{}; j + 1 < n - k; ++j)
                    for (std::size_t i{}; i <= dim; ++i)
                        w[j][i] = (w[j][i] + w[j + 1][i]) / T(2);
            }
            finish_bezier(l);
            finish_bezier(r);
            return {std::move(l), std::move(r)};
        }
    };
}
