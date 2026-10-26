#pragma once
/**
 * @file curvebreaks.h
 * @brief The spans of the parameter range of a curve on which it is smooth: cut at the knots of a B-spline, at the
 * joins of a CurveComposite, and at those of the curves under a CurveTrimmed or a CurveReversed.
 */
#include <gbs/bscurve.h>
#include <gbs/curvecomposite.h>
#include <gbs/curvetrimmed.h>
#include <gbs/curveReversed.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <vector>

namespace gbs
{
    namespace curve_breaks_detail
    {
        /// The parameters in ]a, b[ where c may not be smooth, in any order.
        template <typename T, std::size_t dim>
        void inner_breaks(const Curve<T, dim> &c, T a, T b, std::vector<T> &out)
        {
            auto keep = [&](T u) {
                if (a < u && u < b)
                    out.push_back(u);
            };
            if (const auto *bs = dynamic_cast<const BSCurveGeneral<T, dim, false> *>(&c))
                std::ranges::for_each(bs->knots(), keep);
            else if (const auto *rbs = dynamic_cast<const BSCurveGeneral<T, dim, true> *>(&c))
                std::ranges::for_each(rbs->knots(), keep);
            else if (const auto *cc = dynamic_cast<const CurveComposite<T, dim> *>(&c))
            {
                // a curve covers [u, u + its range] of the composite
                T u = cc->bounds()[0];
                for (const auto &sub : cc->curves())
                {
                    const auto [s0, s1] = sub->bounds();
                    std::vector<T> local;
                    inner_breaks(*sub, s0, s1, local);
                    std::ranges::for_each(local, [&](T s) { keep(u + (s - s0)); });
                    u += s1 - s0;
                    keep(u); // the join (the last one is the end)
                }
            }
            else if (const auto *ct = dynamic_cast<const CurveTrimmed<T, dim> *>(&c))
                inner_breaks(*ct->basisCurve(), a, b, out);
            else if (const auto *cr = dynamic_cast<const CurveReversed<T, dim> *>(&c))
            {
                const auto [u1, u2] = cr->basisCurve()->bounds();
                std::vector<T> basis;
                inner_breaks(*cr->basisCurve(), u1 + u2 - b, u1 + u2 - a, basis);
                std::ranges::for_each(basis, [&](T s) { out.push_back(u1 + u2 - s); });
            }
        }
    } // namespace curve_breaks_detail

    /**
     * @brief The ends of the spans of [a, b] on which a curve is smooth: a, the parameters between where it may not
     * be, in increasing order, then b.
     *
     * A curve may not be smooth at the knots of a B-spline (rational or not), at the joins of a CurveComposite (and
     * at the breaks of its curves), at the breaks of the curve under a CurveTrimmed or a CurveReversed. Any other
     * curve is taken as smooth on [a, b]. The knots outside ]a, b[ are left out, also those a B-spline trimmed
     * without dropping its knots keeps beyond its bounds.
     */
    template <typename T, std::size_t dim>
    auto smooth_spans(const Curve<T, dim> &c, T a, T b) -> std::vector<T>
    {
        std::vector<T> inner;
        curve_breaks_detail::inner_breaks(c, a, b, inner);
        std::ranges::sort(inner);
        const auto [last, end] = std::ranges::unique(inner);
        inner.erase(last, end);
        std::vector<T> spans{a};
        spans.append_range(inner);
        spans.push_back(b);
        return spans;
    }

    /// smooth_spans of a curve over its bounds.
    template <typename T, std::size_t dim>
    auto smooth_spans(const Curve<T, dim> &c) -> std::vector<T>
    {
        const auto [a, b] = c.bounds();
        return smooth_spans(c, a, b);
    }
} // namespace gbs
