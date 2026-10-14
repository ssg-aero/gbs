#pragma once
/**
 * @file box.h
 * @brief Axis-aligned bounding boxes, the filter of the subdivision intersection algorithms.
 */
#include <gbs/gbslib.h>

#include <algorithm>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <limits>
#include <utility>

namespace gbs
{
    /// Axis-aligned box in `dim` dimensions; empty until a point is added.
    template <std::floating_point T, std::size_t dim>
    struct Box
    {
        point<T, dim> min;
        point<T, dim> max;

        Box()
        {
            min.fill(std::numeric_limits<T>::max());
            max.fill(std::numeric_limits<T>::lowest());
        }

        [[nodiscard]] bool empty() const noexcept { return min[0] > max[0]; }

        void add(const point<T, dim> &p) noexcept
        {
            for (std::size_t i{}; i < dim; ++i)
            {
                min[i] = std::min(min[i], p[i]);
                max[i] = std::max(max[i], p[i]);
            }
        }

        void add(const Box &b) noexcept
        {
            if (b.empty())
                return;
            add(b.min);
            add(b.max);
        }

        /// Grows the box by `d` in every direction.
        void inflate(T d) noexcept
        {
            for (std::size_t i{}; i < dim; ++i)
            {
                min[i] -= d;
                max[i] += d;
            }
        }

        /// Do the boxes overlap, or are they closer than `gap`?
        [[nodiscard]] bool overlaps(const Box &b, T gap = T(0)) const noexcept
        {
            if (empty() || b.empty())
                return false;
            for (std::size_t i{}; i < dim; ++i)
                if (b.min[i] > max[i] + gap || min[i] > b.max[i] + gap)
                    return false;
            return true;
        }

        /// Is `b` inside this box grown by `gap`?
        [[nodiscard]] bool contains(const Box &b, T gap = T(0)) const noexcept
        {
            if (empty() || b.empty())
                return false;
            for (std::size_t i{}; i < dim; ++i)
                if (b.min[i] < min[i] - gap || b.max[i] > max[i] + gap)
                    return false;
            return true;
        }

        [[nodiscard]] bool contains(const point<T, dim> &p, T gap = T(0)) const noexcept
        {
            for (std::size_t i{}; i < dim; ++i)
                if (p[i] < min[i] - gap || p[i] > max[i] + gap)
                    return false;
            return true;
        }

        /// Length of the diagonal, 0 for an empty box.
        [[nodiscard]] T diagonal() const noexcept
        {
            if (empty())
                return T(0);
            T d{};
            for (std::size_t i{}; i < dim; ++i)
                d += (max[i] - min[i]) * (max[i] - min[i]);
            return std::sqrt(d);
        }
    };

    /// Does the segment [a, b] meet the box grown by `gap`? (slab test)
    template <std::floating_point T, std::size_t dim>
    [[nodiscard]] bool segment_overlaps(const point<T, dim> &a, const point<T, dim> &b, const Box<T, dim> &box, T gap = T(0))
    {
        if (box.empty())
            return false;
        T lo{0}, hi{1};
        for (std::size_t i{}; i < dim; ++i)
        {
            const T mn = box.min[i] - gap, mx = box.max[i] + gap, d = b[i] - a[i];
            if (d == T(0))
            {
                if (a[i] < mn || a[i] > mx)
                    return false;
                continue;
            }
            T t1 = (mn - a[i]) / d, t2 = (mx - a[i]) / d;
            if (t1 > t2)
                std::swap(t1, t2);
            lo = std::max(lo, t1);
            hi = std::min(hi, t2);
            if (lo > hi)
                return false;
        }
        return true;
    }

    /// Box of a range of points.
    template <std::floating_point T, std::size_t dim, typename Range>
    [[nodiscard]] auto box_of(const Range &points) -> Box<T, dim>
    {
        Box<T, dim> b;
        for (const auto &p : points)
            b.add(p);
        return b;
    }
}
