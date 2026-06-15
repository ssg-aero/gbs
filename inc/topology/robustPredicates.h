#pragma once

// Exact-sign 2D geometric predicates (adaptive precision).
//
// This is a templated C++ port of the orientation predicate from Jonathan
// Richard Shewchuk's public-domain "Routines for Arbitrary Precision
// Floating-Point Arithmetic and Fast Robust Geometric Predicates"
// (https://www.cs.cmu.edu/~quake/robust.html). `orient2d` returns a value whose
// SIGN is the exact sign of the 2x2 determinant
//
//     | ax-cx  ay-cy |
//     | bx-cx  by-cy |
//
// i.e. > 0 iff (a,b,c) turns counter-clockwise, regardless of floating-point
// rounding. A fast floating-point filter handles the common, well-separated
// case; only near-degenerate inputs fall back to the exact expansion path. The
// sign is therefore identical on every platform/compiler, which is what the
// half-edge tessellation needs to build a robust, deterministic mesh
// (cf. issue #48).
//
// The expansion arithmetic relies on IEEE-754 round-to-nearest and on the
// intermediate products/sums NOT being fused (an FMA contraction of `a*b - c*d`
// would defeat the error-free transformations). We disable contraction for this
// translation-unit region; the algorithm is otherwise standard-conforming.

#if defined(__clang__) || defined(__GNUC__)
#pragma STDC FP_CONTRACT OFF
#endif

#include <array>
#include <concepts>
#include <limits>
#include <cmath>

namespace gbs
{
namespace robust_predicates
{
    /// Compile-time constants of Shewchuk's scheme, derived from the type's
    /// machine epsilon. `epsilon` is half an ulp of 1 (round-to-nearest), the
    /// `splitter` is 2^ceil(p/2)+1 used to split a mantissa in two halves, and
    /// the error bounds are the polynomials given in predicates.c.
    template <std::floating_point T>
    struct constants
    {
        static constexpr T epsilon = std::numeric_limits<T>::epsilon() / T{2};

        static constexpr T compute_splitter()
        {
            // 2^( ceil(digits/2) ) + 1, computed exactly in floating point.
            T s{1};
            for (int i = 0; i < (std::numeric_limits<T>::digits + 1) / 2; ++i)
            {
                s *= T{2};
            }
            return s + T{1};
        }

        static constexpr T splitter        = compute_splitter();
        static constexpr T resulterrbound  = (T{3} + T{8} * epsilon) * epsilon;
        static constexpr T ccwerrboundA    = (T{3} + T{16} * epsilon) * epsilon;
        static constexpr T ccwerrboundB    = (T{2} + T{12} * epsilon) * epsilon;
        static constexpr T ccwerrboundC    = (T{9} + T{64} * epsilon) * epsilon * epsilon;
    };

    // --- Error-free transformations (Shewchuk / Dekker / Knuth) --------------

    /// x = fl(a+b), y = a+b-x exactly. No magnitude assumption.
    template <std::floating_point T>
    inline void two_sum(T a, T b, T &x, T &y)
    {
        x = a + b;
        T bvirt  = x - a;
        T avirt  = x - bvirt;
        T bround = b - bvirt;
        T around = a - avirt;
        y = around + bround;
    }

    /// Tail of a-b given the already-rounded difference x = fl(a-b).
    template <std::floating_point T>
    inline void two_diff_tail(T a, T b, T x, T &y)
    {
        T bvirt  = a - x;
        T avirt  = x + bvirt;
        T bround = bvirt - b;
        T around = a - avirt;
        y = around + bround;
    }

    /// x = fl(a-b), y = a-b-x exactly.
    template <std::floating_point T>
    inline void two_diff(T a, T b, T &x, T &y)
    {
        x = a - b;
        two_diff_tail(a, b, x, y);
    }

    /// x = fl(a+b), y tail; requires |a| >= |b|.
    template <std::floating_point T>
    inline void fast_two_sum(T a, T b, T &x, T &y)
    {
        x = a + b;
        T bvirt = x - a;
        y = b - bvirt;
    }

    /// Split a p-bit value into two non-overlapping halves.
    template <std::floating_point T>
    inline void split(T a, T &ahi, T &alo)
    {
        T c    = constants<T>::splitter * a;
        T abig = c - a;
        ahi    = c - abig;
        alo    = a - ahi;
    }

    /// x = fl(a*b), y = a*b-x exactly.
    template <std::floating_point T>
    inline void two_product(T a, T b, T &x, T &y)
    {
        x = a * b;
        T ahi, alo, bhi, blo;
        split(a, ahi, alo);
        split(b, bhi, blo);
        T err1 = x - (ahi * bhi);
        T err2 = err1 - (alo * bhi);
        T err3 = err2 - (ahi * blo);
        y = (alo * blo) - err3;
    }

    /// (a1,a0) - b  ->  (x2,x1,x0)
    template <std::floating_point T>
    inline void two_one_diff(T a1, T a0, T b, T &x2, T &x1, T &x0)
    {
        T i;
        two_diff(a0, b, i, x0);
        two_sum(a1, i, x2, x1);
    }

    /// (a1,a0) - (b1,b0)  ->  (x3,x2,x1,x0)
    template <std::floating_point T>
    inline void two_two_diff(T a1, T a0, T b1, T b0, T &x3, T &x2, T &x1, T &x0)
    {
        T j, z0;
        two_one_diff(a1, a0, b0, j, z0, x0);
        two_one_diff(j, z0, b1, x3, x2, x1);
    }

    /// Sum of the components of an expansion (its rounded value).
    template <std::floating_point T>
    inline T estimate(int elen, const T *e)
    {
        T Q = e[0];
        for (int i = 1; i < elen; ++i)
        {
            Q += e[i];
        }
        return Q;
    }

    /// h = e + f, both non-overlapping increasing-magnitude expansions, with
    /// zero components eliminated. |h| <= elen+flen. The one-past reads of the
    /// canonical C version are guarded so we never index out of bounds.
    template <std::floating_point T>
    inline int fast_expansion_sum_zeroelim(int elen, const T *e, int flen, const T *f, T *h)
    {
        T Q, Qnew, hh;
        int eindex = 0, findex = 0, hindex = 0;
        T enow = e[0];
        T fnow = f[0];

        if ((fnow > enow) == (fnow > -enow))
        {
            Q = enow;
            ++eindex;
            if (eindex < elen) enow = e[eindex];
        }
        else
        {
            Q = fnow;
            ++findex;
            if (findex < flen) fnow = f[findex];
        }

        if ((eindex < elen) && (findex < flen))
        {
            if ((fnow > enow) == (fnow > -enow))
            {
                fast_two_sum(enow, Q, Qnew, hh);
                ++eindex;
                if (eindex < elen) enow = e[eindex];
            }
            else
            {
                fast_two_sum(fnow, Q, Qnew, hh);
                ++findex;
                if (findex < flen) fnow = f[findex];
            }
            Q = Qnew;
            if (hh != T{0}) h[hindex++] = hh;

            while ((eindex < elen) && (findex < flen))
            {
                if ((fnow > enow) == (fnow > -enow))
                {
                    two_sum(Q, enow, Qnew, hh);
                    ++eindex;
                    if (eindex < elen) enow = e[eindex];
                }
                else
                {
                    two_sum(Q, fnow, Qnew, hh);
                    ++findex;
                    if (findex < flen) fnow = f[findex];
                }
                Q = Qnew;
                if (hh != T{0}) h[hindex++] = hh;
            }
        }

        while (eindex < elen)
        {
            two_sum(Q, enow, Qnew, hh);
            ++eindex;
            if (eindex < elen) enow = e[eindex];
            Q = Qnew;
            if (hh != T{0}) h[hindex++] = hh;
        }
        while (findex < flen)
        {
            two_sum(Q, fnow, Qnew, hh);
            ++findex;
            if (findex < flen) fnow = f[findex];
            Q = Qnew;
            if (hh != T{0}) h[hindex++] = hh;
        }

        if ((Q != T{0}) || (hindex == 0))
        {
            h[hindex++] = Q;
        }
        return hindex;
    }

    /// Adaptive stage of orient2d once the cheap filter is inconclusive.
    template <std::floating_point T>
    T orient2d_adapt(const std::array<T, 2> &pa,
                     const std::array<T, 2> &pb,
                     const std::array<T, 2> &pc,
                     T detsum)
    {
        using K = constants<T>;

        T acx = pa[0] - pc[0];
        T bcx = pb[0] - pc[0];
        T acy = pa[1] - pc[1];
        T bcy = pb[1] - pc[1];

        T detleft, detlefttail, detright, detrighttail;
        two_product(acx, bcy, detleft, detlefttail);
        two_product(acy, bcx, detright, detrighttail);

        T B[4];
        two_two_diff(detleft, detlefttail, detright, detrighttail, B[3], B[2], B[1], B[0]);

        T det = estimate(4, B);
        T errbound = K::ccwerrboundB * detsum;
        if ((det >= errbound) || (-det >= errbound))
        {
            return det;
        }

        T acxtail, bcxtail, acytail, bcytail;
        two_diff_tail(pa[0], pc[0], acx, acxtail);
        two_diff_tail(pb[0], pc[0], bcx, bcxtail);
        two_diff_tail(pa[1], pc[1], acy, acytail);
        two_diff_tail(pb[1], pc[1], bcy, bcytail);

        if (acxtail == T{0} && acytail == T{0} && bcxtail == T{0} && bcytail == T{0})
        {
            return det;
        }

        errbound = K::ccwerrboundC * detsum + K::resulterrbound * std::abs(det);
        det += (acx * bcytail + bcy * acxtail) - (acy * bcxtail + bcx * acytail);
        if ((det >= errbound) || (-det >= errbound))
        {
            return det;
        }

        T s1, s0, t1, t0;
        T u[4];

        two_product(acxtail, bcy, s1, s0);
        two_product(acytail, bcx, t1, t0);
        two_two_diff(s1, s0, t1, t0, u[3], u[2], u[1], u[0]);
        T C1[8];
        int C1len = fast_expansion_sum_zeroelim(4, B, 4, u, C1);

        two_product(acx, bcytail, s1, s0);
        two_product(acy, bcxtail, t1, t0);
        two_two_diff(s1, s0, t1, t0, u[3], u[2], u[1], u[0]);
        T C2[12];
        int C2len = fast_expansion_sum_zeroelim(C1len, C1, 4, u, C2);

        two_product(acxtail, bcytail, s1, s0);
        two_product(acytail, bcxtail, t1, t0);
        two_two_diff(s1, s0, t1, t0, u[3], u[2], u[1], u[0]);
        T D[16];
        int Dlen = fast_expansion_sum_zeroelim(C2len, C2, 4, u, D);

        return D[Dlen - 1];
    }
} // namespace robust_predicates

/**
 * @brief Exact-sign orientation of the triangle (a,b,c).
 *
 * Returns a value with the same sign as the exact 2x2 determinant: positive iff
 * (a,b,c) is counter-clockwise, negative iff clockwise, exactly zero iff the
 * three points are collinear. The magnitude approximates the signed area times
 * two but is only guaranteed for its sign. Unlike a raw floating-point
 * determinant, the sign is robust and reproducible across platforms.
 *
 * @tparam T Floating-point type.
 */
template <std::floating_point T>
T orient2d(const std::array<T, 2> &pa, const std::array<T, 2> &pb, const std::array<T, 2> &pc)
{
    using K = robust_predicates::constants<T>;

    T detleft  = (pa[0] - pc[0]) * (pb[1] - pc[1]);
    T detright = (pa[1] - pc[1]) * (pb[0] - pc[0]);
    T det      = detleft - detright;

    T detsum;
    if (detleft > T{0})
    {
        if (detright <= T{0})
            return det;
        detsum = detleft + detright;
    }
    else if (detleft < T{0})
    {
        if (detright >= T{0})
            return det;
        detsum = -detleft - detright;
    }
    else
    {
        return det;
    }

    T errbound = K::ccwerrboundA * detsum;
    if ((det >= errbound) || (-det >= errbound))
    {
        return det;
    }

    return robust_predicates::orient2d_adapt(pa, pb, pc, detsum);
}

} // namespace gbs

// Restore the default contraction setting for the rest of the translation unit
// so including this header does not silently change other code's FP behaviour.
#if defined(__clang__) || defined(__GNUC__)
#pragma STDC FP_CONTRACT DEFAULT
#endif
