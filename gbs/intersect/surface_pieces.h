#pragma once
/**
 * @file surface_pieces.h
 * @brief Patches of a surface for subdivision algorithms: a box enclosing each patch, its
 * flatness, and the split of a patch in two halves along u or v.
 *
 * A B-spline surface (rational or not) is cut into Bézier patches: the box of the Cartesian
 * poles of a patch encloses it, and a patch is halved exactly by de Casteljau on its rows or
 * columns. Any other surface (primitives, revolution, offset…) is sampled on a grid: its box is
 * the box of the samples, grown by a margin estimated from the deviation of the surface at the
 * centre of each cell to the mean of the cell corners. That box is not guaranteed.
 */
#include <gbs/intersect/curve_pieces.h>
#include <gbs/bssurf.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

#ifdef GBS_USE_MODULES
    import vecop;
#endif

namespace gbs
{
    /// A patch of a surface on the parameter box [u1, u2] x [v1, v2].
    template <std::floating_point T>
    struct SurfacePiece
    {
        std::array<T, 4> uv{}; ///< u1, u2, v1, v2
        /// Bézier patch: homogeneous poles (w x, w), u fastest; sampled patch: (x, 1) on a grid
        std::vector<std::array<T, 4>> poles;
        std::size_t nu{}, nv{}; ///< number of poles (or samples) along u and v
        bool bezier{};
        Box<T, 3> box;
        T flatness{}; ///< bound of the distance of the patch to the plane of its corners

        [[nodiscard]] T size() const { return box.diagonal(); }
        [[nodiscard]] std::array<T, 2> center() const { return {(uv[0] + uv[1]) / T(2), (uv[2] + uv[3]) / T(2)}; }
        [[nodiscard]] point<T, 3> pole(std::size_t i, std::size_t j) const
        {
            return intersect_detail::cartesian<T, 3>(poles[i + j * nu]);
        }
    };

    /// Cuts a surface into patches and halves them on demand.
    template <std::floating_point T>
    class SurfacePieces
    {
        const Surface<T, 3> &m_srf;
        std::size_t m_samples;

        // plane of the corners: centroid and normal from the diagonals; distance of the points to it
        static void finish(SurfacePiece<T> &pc, T margin)
        {
            pc.box = Box<T, 3>{};
            point<T, 3> c{};
            for (std::size_t k{}; k < pc.poles.size(); ++k)
            {
                const auto x = intersect_detail::cartesian<T, 3>(pc.poles[k]);
                pc.box.add(x);
                c = c + x;
            }
            c = c / T(pc.poles.size());
            const auto p00 = pc.pole(0, 0), p10 = pc.pole(pc.nu - 1, 0), p01 = pc.pole(0, pc.nv - 1), p11 = pc.pole(pc.nu - 1, pc.nv - 1);
            auto n = cross(p11 - p00, p01 - p10);
            const auto nn = norm(n);
            pc.flatness = T(0);
            for (const auto &h : pc.poles)
            {
                const auto x = intersect_detail::cartesian<T, 3>(h);
                pc.flatness = std::max(pc.flatness, nn > T(0) ? std::abs((x - c) * n) / nn : norm(x - c));
            }
            pc.box.inflate(margin);
            pc.flatness += margin;
        }

        template <typename Bs>
        auto bezier_pieces(Bs s, const std::array<T, 4> &box) -> std::vector<SurfacePiece<T>>
        {
            s.trimU(box[0], box[1]);
            s.trimV(box[2], box[3]);
            const auto p = s.degreeU(), q = s.degreeV();
            const auto ku = s.knotsU(), kv = s.knotsV();
            const auto su = ku.size() - 1, sv = kv.size() - 1;
            std::vector<SurfacePiece<T>> pieces;
            // Bézier decomposition row by row along u, then column by column along v (The NURBS
            // Book, A5.7 by A5.6 on each row and column), every knot raised at once
            using Pole = std::remove_cvref_t<decltype(s.poles()[0])>;
            const auto n0U = s.nPolesU(), n0V = s.nPolesV();
            const auto nU = su * p + 1, nV = sv * q + 1;
            std::vector<Pole> rows(nU * n0V);
            {
                std::vector<Pole> row(n0U);
                for (std::size_t j = 0; j < n0V; ++j)
                {
                    for (std::size_t i = 0; i < n0U; ++i)
                        row[i] = s.poles()[i + j * n0U];
                    const auto r = bezier_poles(s.knotsFlatsU(), row, p);
                    if (!r || r->size() != nU)
                        return pieces; // not clamped, or discontinuous: the caller samples the surface
                    for (std::size_t i = 0; i < nU; ++i)
                        rows[i + j * nU] = (*r)[i];
                }
            }
            std::vector<Pole> poles(nU * nV);
            {
                std::vector<Pole> col(n0V);
                for (std::size_t i = 0; i < nU; ++i)
                {
                    for (std::size_t j = 0; j < n0V; ++j)
                        col[j] = rows[i + j * nU];
                    const auto r = bezier_poles(s.knotsFlatsV(), col, q);
                    if (!r || r->size() != nV)
                        return pieces;
                    for (std::size_t j = 0; j < nV; ++j)
                        poles[i + j * nU] = (*r)[j];
                }
            }
            for (std::size_t b = 0; b < sv; ++b)
                for (std::size_t a = 0; a < su; ++a)
                {
                    SurfacePiece<T> pc{.uv = {ku[a], ku[a + 1], kv[b], kv[b + 1]}, .nu = p + 1, .nv = q + 1, .bezier = true};
                    pc.poles.resize((p + 1) * (q + 1));
                    for (std::size_t j = 0; j <= q; ++j)
                        for (std::size_t i = 0; i <= p; ++i)
                        {
                            const auto &x = poles[(a * p + i) + (b * q + j) * nU];
                            auto &h = pc.poles[i + j * (p + 1)];
                            if constexpr (std::tuple_size_v<std::remove_cvref_t<decltype(x)>> == 4)
                                h = x;
                            else
                                h = {x[0], x[1], x[2], T(1)};
                        }
                    finish(pc, T(0));
                    pieces.push_back(std::move(pc));
                }
            return pieces;
        }

        auto sampled_piece(const std::array<T, 4> &uv) const -> SurfacePiece<T>
        {
            const auto n = m_samples;
            SurfacePiece<T> pc{.uv = uv, .nu = n + 1, .nv = n + 1, .bezier = false};
            pc.poles.resize((n + 1) * (n + 1));
            auto u_at = [&](T k) { return uv[0] + (uv[1] - uv[0]) * k / T(n); };
            auto v_at = [&](T k) { return uv[2] + (uv[3] - uv[2]) * k / T(n); };
            for (std::size_t j = 0; j <= n; ++j)
                for (std::size_t i = 0; i <= n; ++i)
                {
                    const auto x = m_srf.value(i == n ? uv[1] : u_at(T(i)), j == n ? uv[3] : v_at(T(j)));
                    pc.poles[i + j * (n + 1)] = {x[0], x[1], x[2], T(1)};
                }
            // margin: twice the largest deviation of a cell centre to the mean of its corners
            T margin{};
            for (std::size_t j = 0; j < n; ++j)
                for (std::size_t i = 0; i < n; ++i)
                {
                    const auto mid = m_srf.value(u_at(T(i) + T(0.5)), v_at(T(j) + T(0.5)));
                    const auto mean = (pc.pole(i, j) + pc.pole(i + 1, j) + pc.pole(i, j + 1) + pc.pole(i + 1, j + 1)) / T(4);
                    margin = std::max(margin, norm(mid - mean));
                }
            finish(pc, T(2) * margin);
            return pc;
        }

        // de Casteljau at 1/2 on the lines of a grid: along u (rows) or v (columns)
        static auto halve(const SurfacePiece<T> &pc, bool along_u) -> std::pair<SurfacePiece<T>, SurfacePiece<T>>
        {
            SurfacePiece<T> l = pc, r = pc;
            const auto n = along_u ? pc.nu : pc.nv, lines = along_u ? pc.nv : pc.nu;
            auto at = [&](std::size_t line, std::size_t k) { return along_u ? k + line * pc.nu : line + k * pc.nu; };
            std::vector<std::array<T, 4>> w(n);
            for (std::size_t line = 0; line < lines; ++line)
            {
                for (std::size_t k = 0; k < n; ++k)
                    w[k] = pc.poles[at(line, k)];
                for (std::size_t k = 0; k < n; ++k)
                {
                    l.poles[at(line, k)] = w.front();
                    r.poles[at(line, n - 1 - k)] = w[n - 1 - k];
                    for (std::size_t j = 0; j + 1 < n - k; ++j)
                        for (std::size_t i = 0; i < 4; ++i)
                            w[j][i] = (w[j][i] + w[j + 1][i]) / T(2);
                }
            }
            const auto c = pc.center();
            if (along_u)
                l.uv[1] = r.uv[0] = c[0];
            else
                l.uv[3] = r.uv[2] = c[1];
            finish(l, T(0));
            finish(r, T(0));
            return {std::move(l), std::move(r)};
        }

    public:
        /// @param samples number of cells sampled along each direction of a patch of a surface which is not a B-spline
        explicit SurfacePieces(const Surface<T, 3> &srf, std::size_t samples = 6) : m_srf{srf}, m_samples{std::max<std::size_t>(samples, 2)} {}

        [[nodiscard]] auto surface() const -> const Surface<T, 3> & { return m_srf; }

        /// The initial patches of the surface on the box {u1, u2, v1, v2}: its Bézier patches, or one sampled patch.
        [[nodiscard]] auto pieces(const std::array<T, 4> &box) -> std::vector<SurfacePiece<T>>
        {
            std::vector<SurfacePiece<T>> r;
            if (const auto *bs = dynamic_cast<const BSSurface<T, 3> *>(&m_srf))
                r = bezier_pieces(*bs, box);
            else if (const auto *bsr = dynamic_cast<const BSSurfaceRational<T, 3> *>(&m_srf))
                r = bezier_pieces(*bsr, box);
            if (r.empty())
                r.push_back(sampled_piece(box));
            return r;
        }

        /// The two halves of a patch along u or v, at the middle of its parameter range.
        [[nodiscard]] auto split(const SurfacePiece<T> &pc, bool along_u) const -> std::pair<SurfacePiece<T>, SurfacePiece<T>>
        {
            if (pc.bezier)
                return halve(pc, along_u);
            const auto c = pc.center();
            auto a = pc.uv, b = pc.uv;
            if (along_u)
                a[1] = b[0] = c[0];
            else
                a[3] = b[2] = c[1];
            return {sampled_piece(a), sampled_piece(b)};
        }

        /// The direction to split a patch: the one along which it is the longest in 3D, measured
        /// along the lines of poles (end points coincide on a closed patch).
        [[nodiscard]] static bool longer_along_u(const SurfacePiece<T> &pc)
        {
            T lu{}, lv{};
            for (std::size_t j : {std::size_t{0}, pc.nv / 2, pc.nv - 1})
            {
                T l{};
                for (std::size_t i = 0; i + 1 < pc.nu; ++i)
                    l += norm(pc.pole(i + 1, j) - pc.pole(i, j));
                lu = std::max(lu, l);
            }
            for (std::size_t i : {std::size_t{0}, pc.nu / 2, pc.nu - 1})
            {
                T l{};
                for (std::size_t j = 0; j + 1 < pc.nv; ++j)
                    l += norm(pc.pole(i, j + 1) - pc.pole(i, j));
                lv = std::max(lv, l);
            }
            return lu >= lv;
        }
    };
}
