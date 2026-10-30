#pragma once
/**
 * @file polynomial.h
 * @brief Real roots of a polynomial of low degree on an interval.
 */
#include <algorithm>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <vector>

namespace gbs
{
    /// Value of c[0] + c[1] t + ... + c[n] t^n (Horner).
    template <std::floating_point T>
    [[nodiscard]] T polynomial_value(const std::vector<T> &c, T t)
    {
        T r{};
        for (auto it = c.rbegin(); it != c.rend(); ++it)
            r = r * t + *it;
        return r;
    }

    template <std::floating_point T>
    struct PolynomialRoots
    {
        std::vector<T> roots;    ///< sign changes, sorted
        std::vector<T> extremes; ///< local extremes, sorted: double roots are among them
    };

    /**
     * @brief Roots of a polynomial in [a, b], isolated between the roots of its derivative and
     * refined by bisection. Double roots do not change sign: they are local extremes, returned
     * apart for the caller to test. A root at an end of [a, b] may be missed by the rounding:
     * the caller tests the ends.
     */
    template <std::floating_point T>
    [[nodiscard]] auto polynomial_roots(std::vector<T> c, T a, T b) -> PolynomialRoots<T>
    {
        // drop the leading coefficients whose term vanishes on [a, b]: |c_k| h^k against the
        // largest term max |c_i| h^i, h the largest |t| (the coefficients alone mix powers)
        T h = std::max(std::abs(a), std::abs(b));
        if (!(h > T(0)))
            h = T(1);
        std::vector<T> term(c.size());
        T scale{}, hk{1};
        for (std::size_t k{}; k < c.size(); ++k, hk *= h)
        {
            term[k] = std::abs(c[k]) * hk;
            scale = std::max(scale, term[k]);
        }
        while (c.size() > 1 && term[c.size() - 1] <= T(1e-14) * scale)
            c.pop_back();
        PolynomialRoots<T> r;
        if (c.size() < 2)
            return r;
        if (c.size() == 2)
        {
            const T t = -c[0] / c[1];
            if (t >= a && t <= b)
                r.roots.push_back(t);
            return r;
        }
        std::vector<T> d(c.size() - 1);
        for (std::size_t i = 1; i < c.size(); ++i)
            d[i - 1] = T(i) * c[i];
        r.extremes = polynomial_roots(d, a, b).roots;
        std::vector<T> knots{a};
        knots.insert(knots.end(), r.extremes.begin(), r.extremes.end());
        knots.push_back(b);
        for (std::size_t k = 0; k + 1 < knots.size(); ++k)
        {
            T lo = knots[k], hi = knots[k + 1];
            T flo = polynomial_value(c, lo), fhi = polynomial_value(c, hi);
            if (flo == T(0))
            {
                if (r.roots.empty() || r.roots.back() != lo)
                    r.roots.push_back(lo);
                continue;
            }
            if (fhi == T(0))
            {
                if (k + 2 == knots.size())
                    r.roots.push_back(hi);
                continue;
            }
            if ((flo < T(0)) == (fhi < T(0)))
                continue;
            for (int it{}; it < 200 && hi - lo > T(0); ++it)
            {
                const T m = (lo + hi) / T(2);
                if (m <= lo || m >= hi)
                    break;
                const T fm = polynomial_value(c, m);
                if ((fm < T(0)) == (flo < T(0)))
                {
                    lo = m;
                    flo = fm;
                }
                else
                    hi = m;
            }
            r.roots.push_back((lo + hi) / T(2));
        }
        return r;
    }
}
