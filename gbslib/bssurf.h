#pragma once
#include <gbslib/basisfunctions.h>
#include <gbslib/knotsfunctions.h>
#include <gbslib/bscurve.h>
#include <vector>
#include <array>
#include <algorithm>
namespace gbs
{
    /**
     * @brief Géneral BSpline surface class, any kind of precision, space dimension with rational definition capability
     * 
     * @tparam T 
     * @tparam dim 
     */
    template <typename T, size_t dim, bool rational>
    class BSSurfaceGeneral
    {
        size_t m_degU,m_degV;
        std::vector< std::array<T, dim + rational> > m_poles;
        std::vector<T> m_knotsFlatsU;
        std::vector<T> m_knotsFlatsV;

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
                size_t degV
                ) : m_poles(poles),
                              m_knotsFlatsU(knots_flatsU),
                              m_knotsFlatsV(knots_flatsV),
                              m_degU(degU),
                              m_degV(degV)//,
                            //   m_rational(false)

        {
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
                                 m_degV(degV)//,
                                //  m_rational(true)

        {
        }
        /**
         * @brief  Non rational curve evaluation
         * 
         * @param u  : u parameter on surface
         * @param v  : v parameter on surface
         * @param du : u derivative order
         * @param dv : v derivative order
         * @return std::array<T, dim> const 
         */
        virtual auto value(T u, T v, size_t du = 0, size_t dv = 0) const -> std::array<T, dim> = 0;

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

        // auto insertKnotU(T u, size_t m = 1) -> void //Fail saife, i.e. if fails, curve stays in precedent state
        // {
        //     for (auto i = 0; i < m; i++)
        //         insert_knot(u, m_deg, m_knotsFlatsU, m_poles);
        // }

        // auto removeKnot(T u, T tol, size_t m = 1) -> void //Fail saife, i.e. if fails, curve stays in precedent state
        // {
        //     for (auto i = 0; i < m; i++)
        //         remove_knot(u, m_deg, m_knotsFlats, m_poles, tol);
        // }

        auto poles() const noexcept -> const std::vector<std::array<T, dim + rational>> &
        {
            return m_poles;
        }

        auto nPolesU() const noexcept -> size_t
        {
            return m_knotsFlatsU.size() - m_degU - 1;
        }

        auto nPolesV() const noexcept -> size_t
        {
            return m_knotsFlatsV.size() - m_degV - 1;
        }
    };

    template <typename T, size_t dim>
    class BSSurface : public BSSurfaceGeneral<T, dim, false>
    {
    public:
        BSSurface(const std::vector<std::array<T, dim>> &poles,
                const std::vector<T> &knots_flatsU,
                const std::vector<T> &knots_flatsV,
                size_t degU,
                size_t degV
                ) : BSSurfaceGeneral<T, dim, false>(poles, knots_flatsU,knots_flatsV, degU,degV) {}
        virtual auto value(T u, T v, size_t du = 0, size_t dv = 0) const -> std::array<T, dim> override
        {
            return gbs::eval_value_simple(u, v, knotsFlatsU(), knotsFlatsV(), poles(), degreeU(), degreeV(), du, dv);
        }
    };

    template <typename T, size_t dim>
    class BSSurfaceRational : public BSSurfaceGeneral<T, dim, true>
    {
    public:
        BSSurfaceRational(  const std::vector<std::array<T, dim + 1>> &poles,
                            const std::vector<T> &knots_flatsU,
                            const std::vector<T> &knots_flatsV,
                            size_t degU,
                            size_t degV
                            ) : BSSurfaceGeneral<T, dim, true>(poles, knots_flatsU,knots_flatsV, degU,degV) {}
        virtual auto value(T u, T v, size_t du = 0, size_t dv = 0) const -> std::array<T, dim> override
        {
            if (du == 0 && dv == 0)
            {
                return weight_projection( gbs::eval_value_simple(u, v, knotsFlatsU(), knotsFlatsV(), poles(), degreeU(), degreeV(), du, dv));
            }
            else
            {
                throw std::exception("not implemented");
            }
        }

        auto polesProjected() const -> points_vector<T,dim>
        {
            points_vector<T,dim> poles_(poles().size());
            std::transform(poles().begin(),poles().end(),poles_.begin(),
            [](const auto &p){return weight_projection(p);});
            return poles_;
        }

        auto weights() const -> std::vector<T>
        {
            std::vector<T> weights_(poles().size());
            std::transform(poles().begin(),poles().end(),weights_.begin(),
            [](const auto &p){return p.back();});
            return weights_;
        }
    };

} // namespace gbs