#pragma once

/**
 * @file placement.h
 * @brief Native BREP core: rigid placement (rotation and translation).
 *
 * Used by copy_shape (gbs-brep/edit.h) and by the instance provenance of the
 * model (gbs-brep/model.h). Architecture notes docs/sources/design/step_pr09_assemblies.md
 * and step_p2_instances.md.
 */

#include <gbs/gbslib.h>

#include <array>
#include <cmath>

namespace gbs::brep
{
    /**
     * @brief Rigid placement p -> R p + t (rotation R, translation t), as given by a STEP
     * axis2_placement_3d or an assembly transformation.
     */
    template <std::floating_point T>
    struct Placement
    {
        std::array<std::array<T, 3>, 3> r{{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}}; ///< rows of R
        point<T, 3> t{};

        /// Placement of the frame (origin, x, y, z): local coordinates -> global ones.
        [[nodiscard]] static Placement frame(const point<T, 3> &o, const point<T, 3> &x, const point<T, 3> &y, const point<T, 3> &z)
        {
            Placement p;
            for (std::size_t i = 0; i < 3; ++i)
                p.r[i] = {x[i], y[i], z[i]};
            p.t = o;
            return p;
        }

        [[nodiscard]] point<T, 3> apply(const point<T, 3> &p) const
        {
            point<T, 3> q;
            for (std::size_t i = 0; i < 3; ++i)
                q[i] = r[i][0] * p[0] + r[i][1] * p[1] + r[i][2] * p[2] + t[i];
            return q;
        }

        /// this o b: b applied first
        [[nodiscard]] Placement operator*(const Placement &b) const
        {
            Placement c;
            for (std::size_t i = 0; i < 3; ++i)
                for (std::size_t j = 0; j < 3; ++j)
                    c.r[i][j] = r[i][0] * b.r[0][j] + r[i][1] * b.r[1][j] + r[i][2] * b.r[2][j];
            c.t = apply(b.t);
            return c;
        }

        /// Inverse of a rigid placement (R orthonormal).
        [[nodiscard]] Placement inverse() const
        {
            Placement c;
            for (std::size_t i = 0; i < 3; ++i)
                for (std::size_t j = 0; j < 3; ++j)
                    c.r[i][j] = r[j][i];
            for (std::size_t i = 0; i < 3; ++i)
                c.t[i] = -(c.r[i][0] * t[0] + c.r[i][1] * t[1] + c.r[i][2] * t[2]);
            return c;
        }

        [[nodiscard]] bool is_identity(T eps = T(1e-14)) const
        {
            for (std::size_t i = 0; i < 3; ++i)
            {
                if (std::abs(t[i]) > eps)
                    return false;
                for (std::size_t j = 0; j < 3; ++j)
                    if (std::abs(r[i][j] - T(i == j)) > eps)
                        return false;
            }
            return true;
        }
    };

}
