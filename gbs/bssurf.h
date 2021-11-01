#pragma once
#include <gbs/basisfunctions.h>
#include <gbs/knotsfunctions.h>
#include <gbs/bscurve.h>
#include <vector>
#include <array>
#include <algorithm>
namespace gbs
{
    /**
     * @brief 
     * 
     * @tparam T 
     * @tparam dim 
     * @param poles 
     * @param n_new_poles_u : number of new u poles
     * @return points_vector<T, dim> 
     */
    template <typename T, size_t dim>
    auto inverted_uv_poles(const points_vector<T, dim> &poles, size_t n_new_poles_u) -> points_vector<T, dim>
    {
        auto n_poles_v = poles.size() / n_new_poles_u;
        points_vector<T, dim> poles_t(poles.size());
        for (int i = 0; i < n_new_poles_u; i++)
        {
            for (int j = 0; j < n_poles_v; j++)
            {
                poles_t[i + n_new_poles_u * j] = poles[j + n_poles_v * i];
            }
        }
        return poles_t;
    }
    /**
     * @brief 
     * 
     * @tparam T 
     * @tparam dim 
     * @param poles 
     * @param n_new_poles_u : number of new u poles
     */
    template <typename T, size_t dim>
    auto invert_uv_poles(points_vector<T, dim> &poles, size_t n_new_poles_u) -> void
    {
        poles = std::move(inverted_uv_poles(poles, n_new_poles_u));
    }

    template <typename T, size_t dim>
    class Surface
    {
    public:
        /**
         * @brief  Surface evaluation at parameters {u,v}
         * 
         * @param u  : u parameter on surface
         * @param v  : v parameter on surface
         * @param du : u derivative order
         * @param dv : v derivative order
         * @return point<T, dim> const 
         */
        virtual auto value(T u, T v, size_t du = 0, size_t dv = 0) const -> point<T, dim> = 0;
        /**
         * @brief return surface's bounds {U1,U2,V1,V2}
         * 
         * @return point<T, dim>
         */
        virtual auto bounds() const -> std::array<T, 4> = 0;

        auto operator()(T u, T v, size_t du = 0, size_t dv = 0) const -> point<T, dim> { return value(u, v, du, dv); };
    };

    /**
     * @brief Géneral BSpline surface class, any kind of precision, space dimension with rational definition capability
     * 
     * @tparam T 
     * @tparam dim 
     */
    template <typename T, size_t dim, bool rational>
    class BSSurfaceGeneral : public Surface<T, dim>
    {
        size_t m_degU, m_degV;
        std::vector<std::array<T, dim + rational>> m_poles;
        std::vector<T> m_knotsFlatsU;
        std::vector<T> m_knotsFlatsV;
        std::array<T, 4> m_bounds;

    public:
        /**
     * @brief Construct a new BSSurface object, non rational definition
     * 
     * @param poles 
     * @param knots_flatsU 
     * @param knots_flatsV 
     * @param degU 
     * @param degV 
     */
        BSSurfaceGeneral(const std::vector<std::array<T, dim + rational>> &poles,
                         const std::vector<T> &knots_flatsU,
                         const std::vector<T> &knots_flatsV,
                         size_t degU,
                         size_t degV) : m_poles(poles),
                                        m_knotsFlatsU(knots_flatsU),
                                        m_knotsFlatsV(knots_flatsV),
                                        m_degU(degU),
                                        m_degV(degV),
                                        m_bounds{knots_flatsU.front(), knots_flatsU.back(), knots_flatsV.front(), knots_flatsV.back()}
        //,
        //   m_rational(false)

        {
            if (nPolesU()*nPolesV()!=m_poles.size())
                throw std::exception("BSpline Surface constructor error.");
            if (!check_curve(nPolesU(), m_knotsFlatsU, degU))
                throw std::exception("BSpline Surface constructor error.");
            if (!check_curve(nPolesV(), m_knotsFlatsV, degV))
                throw std::exception("BSpline Surface constructor error.");
        }
        /**
     * @brief Construct a new BSSurface object, rational definition
     * 
     * @param poles 
     * @param weights 
     * @param knots_flatsU 
     * @param knots_flatsV 
     * @param degU 
     * @param degV 
     */
        BSSurfaceGeneral(const std::vector<std::array<T, dim>> &poles,
                         const std::vector<T> &weights,
                         const std::vector<T> &knots_flatsU,
                         const std::vector<T> &knots_flatsV,
                         size_t degU,
                         size_t degV) : m_poles(merge_weights(poles, weights)),
                                        m_knotsFlatsU(knots_flatsU),
                                        m_knotsFlatsV(knots_flatsV),
                                        m_degU(degU),
                                        m_degV(degV),
                                        m_bounds{knots_flatsU.front(), knots_flatsU.back(), knots_flatsV.front(), knots_flatsV.back()}
        //,
        //  m_rational(true)

        {
            if (nPolesU()*nPolesV()!=m_poles.size())
                throw std::exception("BSpline Surface constructor error.");
            if (!check_curve(nPolesU(), m_knotsFlatsU, degU))
                throw std::exception("BSpline Surface constructor error.");
            if (!check_curve(nPolesV(), m_knotsFlatsV, degV))
                throw std::exception("BSpline Surface constructor error.");
        }

        BSSurfaceGeneral() = default;
        BSSurfaceGeneral(const BSSurfaceGeneral<T, dim, rational> &) = default;
        // BSSurfaceGeneral<T, dim, rational> &operator=(BSSurfaceGeneral<T, dim, rational> &srf) const = default;

        auto degreeU() const noexcept -> size_t
        {
            return m_degU;
        }

        auto degreeV() const noexcept -> size_t
        {
            return m_degV;
        }

        auto knotsFlatsU() const noexcept -> const std::vector<T> &
        {
            return m_knotsFlatsU;
        }

        auto knotsFlatsV() const noexcept -> const std::vector<T> &
        {
            return m_knotsFlatsV;
        }

        auto insertKnotU(T u, size_t m = 1) -> void //Fail safe, i.e. if fails, surf stays in previous state
        {
            auto nV = nPolesV();

            points_vector<T, dim + rational> poles_;
            std::vector<T> knots_flatsU_;

            for (auto j = 0; j < nV; j++)
            {
                auto polesU_ = polesU(j);
                knots_flatsU_ = m_knotsFlatsU;

                for (auto i = 0; i < m; i++)
                {
                    insert_knot(u, m_degU, knots_flatsU_, polesU_);
                }

                poles_.insert(poles_.end(), polesU_.begin(), polesU_.end());
            }

            m_knotsFlatsU = std::move(knots_flatsU_);
            m_poles = std::move(poles_);
        }

        auto insertKnotV(T u, size_t m = 1) -> void //Fail safe, i.e. if fails, surf stays in previous state
        {
            auto nU = nPolesU();

            points_vector<T, dim + rational> poles_;
            std::vector<T> knots_flatsV_;

            for (auto j = 0; j < nU; j++)
            {
                auto polesV_ = polesV(j);
                knots_flatsV_ = m_knotsFlatsV;

                for (auto i = 0; i < m; i++)
                {
                    insert_knot(u, m_degV, knots_flatsV_, polesV_);
                }

                poles_.insert(poles_.end(), polesV_.begin(), polesV_.end());
            }

            invert_uv_poles(poles_, nPolesU());

            m_knotsFlatsV = std::move(knots_flatsV_);
            m_poles = std::move(poles_);
        }

        // auto removeKnot(T u, T tol, size_t m = 1) -> void //Fail safe, i.e. if fails, curve stays in previous state
        // {
        //     for (auto i = 0; i < m; i++)
        //         remove_knot(u, m_deg, m_knotsFlats, m_poles, tol);
        // }

        auto poles() const noexcept -> const points_vector<T, dim + rational> &
        {
            return m_poles;
        }
        /*
            auto polesV(size_t j) const noexcept -> const points_vector<T, dim + rational>
            {
                auto nV = nPolesV();
                points_vector<T, dim + rational> poles_(nV);
                auto beg_ = std::next(m_poles.begin(), j * nV);
                auto end_ = std::next(m_poles.begin(), (j + 1) * nV);
                std::copy(std::execution::par, beg_, end_, poles_.begin());

                return poles_;
            }

            auto polesU(size_t j) const noexcept -> const points_vector<T, dim + rational>
            {
                auto nU = nPolesU();
                auto nV = nPolesV();
                points_vector<T, dim + rational> poles_(nU);
                auto it = std::next(m_poles.begin(), j );
                for(auto i = 0;  i < nU; i++ )
                {
                    poles_[i] = *it;
                    it = std::next(it,nV);
                }

                return poles_;
            }
*/
        auto polesU(size_t j) const noexcept -> const points_vector<T, dim + rational>
        {
            auto nU = nPolesU();
            auto nV = nPolesV();
            points_vector<T, dim + rational> poles_(nU);
            auto beg_ = std::next(m_poles.begin(), j * nU);
            auto end_ = (j < nV) ? std::next(m_poles.begin(), (j + 1) * nU) : m_poles.end();
            std::copy(std::execution::par, beg_, end_, poles_.begin());

            return poles_;
        }

        auto polesV(size_t j) const noexcept -> const points_vector<T, dim + rational>
        {
            auto nU = nPolesU();
            auto nV = nPolesV();
            points_vector<T, dim + rational> poles_(nV);
            for (auto i = 0; i < nV; i++)
            {
                poles_[i] = m_poles[j + i * nU];
            }

            return poles_;
        }

        auto nPolesU() const noexcept -> size_t
        {
            return m_knotsFlatsU.size() - m_degU - 1;
        }

        auto nPolesV() const noexcept -> size_t
        {
            return m_knotsFlatsV.size() - m_degV - 1;
        }
        constexpr auto poles_begin() const noexcept { return m_poles.begin(); }
        constexpr auto poles_end() const noexcept { return m_poles.end(); }
        /**
         * @brief Copy Poles, throw std::length_error is thrown if lengths are not the same
         * 
         * @param poles 
         */
        auto copyPoles(const std::vector<std::array<T, dim + rational>> &poles) -> void
        {
            if (poles.size() != m_poles.size())
            {
                throw std::length_error("BSSurfaceGeneral: wrong pole vector length.");
            }
            m_poles = poles;
        }
        /**
         * @brief Move pole vector, , throw std::length_error is thrown if lengths are not the same
         * 
         * @param poles 
         */
        auto movePoles(std::vector<std::array<T, dim + rational>> &poles) -> void
        {
            if (poles.size() != m_poles.size())
            {
                throw std::length_error("BSSurfaceGeneral: wrong pole vector length.");
            }
            m_poles = std::move(poles);
        }

        virtual auto bounds() const -> std::array<T, 4> override
        {
            // return {m_knotsFlatsU.front(), m_knotsFlatsU.back(), m_knotsFlatsV.front(), m_knotsFlatsV.back()};
            return m_bounds;
        }
        auto changeUBounds(T k1, T k2) -> void
        {
            gbs::change_bounds(k1, k2, m_knotsFlatsU);
            m_bounds[0] = k1;
            m_bounds[1] = k2;
        }
        auto changeVBounds(T k1, T k2) -> void
        {
            gbs::change_bounds(k1, k2, m_knotsFlatsV);
            m_bounds[2] = k1;
            m_bounds[3] = k2;
        }
        // virtual auto isoU(T u) const -> BSCurveGeneral<T,dim,rational> = 0;

        // virtual auto isoV(T v) const -> BSCurveGeneral<T,dim,rational> = 0;

        auto increaseDegreeU() -> void
        {
            points_vector<T, dim> poles_new;
            std::vector<T> ku_new;
            std::vector<T> ku_old{m_knotsFlatsU};
            auto nv = nPolesV();
            auto nu = nPolesU();
            for (size_t i{}; i < nv; i++)
            {

                points_vector<T, dim> poles{std::next(m_poles.begin(), i * nu), std::next(m_poles.begin(), i * nu + nu)};
                if (i != 0)
                {
                    m_knotsFlatsU = {ku_old};
                }
                gbs::increase_degree(m_knotsFlatsU, poles, m_degU);
                if (i == 0)
                {
                    ku_new = std::move(m_knotsFlatsU);
                }
                poles_new.insert(poles_new.end(), poles.begin(), poles.end());
            }
            m_poles = std::move(poles_new);
            m_knotsFlatsU = std::move(ku_new);
            m_degU++;
        }

        auto increaseDegreeV() -> void
        {
            invertUV();
            increaseDegreeU();
            invertUV();
        }

        auto invertUV() -> void
        {
            invert_uv_poles(m_poles,nPolesV());
            std::swap(m_knotsFlatsU,m_knotsFlatsV);
            std::swap(m_degU,m_degV);
            m_bounds = {m_knotsFlatsU.front(), m_knotsFlatsU.back(), m_knotsFlatsV.front(), m_knotsFlatsV.back()};
        }
        auto reverseU() -> void
        {
            auto k1 = m_knotsFlatsU.front();
            auto k2 = m_knotsFlatsU.back();
            std::reverse(m_knotsFlatsU.begin(), m_knotsFlatsU.end());
            std::transform(m_knotsFlatsU.begin(), m_knotsFlatsU.end(),
                           m_knotsFlatsU.begin(),
                           [&](const auto k_) {
                               return k1 + k2 - k_;
                           });
                        auto nv = nPolesV();
            auto nu = nPolesU();
            for (size_t i{}; i < nv; i++)
            {
                std::reverse(std::next(m_poles.begin(), i * nu), std::next(m_poles.begin(), i * nu + nu));
            }
        }

        auto reverseV() -> void
        {
            invertUV();
            reverseU();
            invertUV();
        }

    };

    template <typename T, size_t dim>
    class BSSurface : public BSSurfaceGeneral<T, dim, false>
    {
    public:
        using BSSurfaceGeneral<T, dim, false>::BSSurfaceGeneral;
        virtual auto value(T u, T v, size_t du = 0, size_t dv = 0) const -> std::array<T, dim> override
        {
            if (u < this->bounds()[0] - knot_eps || u > this->bounds()[1] + knot_eps)
                throw std::exception("BSpline Curve eval out of U bounds error.");
            if (v < this->bounds()[2] - knot_eps || v > this->bounds()[3] + knot_eps)
                throw std::exception("BSpline Curve eval out of V bounds error.");

            return gbs::eval_value_simple(u, v, this->knotsFlatsU(), this->knotsFlatsV(), this->poles(), this->degreeU(), this->degreeV(), du, dv);
        }

        // virtual auto isoU(T u) const override -> BSCurveGeneral<T,dim,rational>
        // {
        //     auto cpy_ {*this};
        //     cpy_.insertKnotU(u,degreeU()-1);

        // }
    };

    template <typename T, size_t dim>
    class BSSurfaceRational : public BSSurfaceGeneral<T, dim, true>
    {
        using BSSurfaceGeneral<T, dim, true>::BSSurfaceGeneral;
    public:
        BSSurfaceRational(const std::vector<std::array<T, dim>> &poles,
                          const std::vector<T> &knots_flatsU,
                          const std::vector<T> &knots_flatsV,
                          size_t degU,
                          size_t degV) : BSSurfaceGeneral<T, dim, true>(add_weights_coord(poles), knots_flatsU, knots_flatsV, degU, degV) {}
        virtual auto value(T u, T v, size_t du = 0, size_t dv = 0) const -> std::array<T, dim> override
        {
            if (u < this->bounds()[0] - knot_eps || u > this->bounds()[1] + knot_eps)
                throw std::exception("BSpline Curve eval out of U bounds error.");
            if (v < this->bounds()[2] - knot_eps || v > this->bounds()[3] + knot_eps)
                throw std::exception("BSpline Curve eval out of V bounds error.");

            if (du == 0 && dv == 0)
            {
                return weight_projection(gbs::eval_value_simple(u, v, this->knotsFlatsU(), this->knotsFlatsV(), this->poles(), this->degreeU(), this->degreeV(), du, dv));
            }
            else
            {
                throw std::exception("not implemented");
            }
        }

        auto polesProjected() const -> points_vector<T, dim>
        {
            points_vector<T, dim> poles_(this->poles().size());
            std::transform(this->poles().begin(), this->poles().end(), poles_.begin(),
                           [](const auto &p)
                           { return weight_projection(p); });
            return poles_;
        }

        auto weights() const -> std::vector<T>
        {
            std::vector<T> weights_(this->poles().size());
            std::transform(this->poles().begin(), this->poles().end(), weights_.begin(),
                           [](const auto &p)
                           { return p.back(); });
            return weights_;
        }
    };

    template <typename T>
    class BSSfunction
    {
        BSSurface<T, 1> srf_;

    public:
        BSSfunction(const BSSurface<T, 1> &crv) : srf_{crv}
        {
        }
        auto operator()(T u, T v, size_t du = 0, size_t dv = 0) const -> T
        {
            return srf_.value(u, v, du, dv)[0];
        }
        auto bounds() -> std::array<T, 4> { return srf_.bounds(); }
    };

    using BSSurface2d_f = BSSurface<float, 2>;
    using BSSurface3d_f = BSSurface<float, 3>;
    using BSSurface2d_d = BSSurface<double, 2>;
    using BSSurface3d_d = BSSurface<double, 3>;
    using BSSurfaceRational2d_f = BSSurfaceRational<float, 2>;
    using BSSurfaceRational3d_f = BSSurfaceRational<float, 3>;
    using BSSurfaceRational2d_d = BSSurfaceRational<double, 2>;
    using BSSurfaceRational3d_d = BSSurfaceRational<double, 3>;

} // namespace gbs