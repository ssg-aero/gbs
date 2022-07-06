#pragma once
#include <gbs/curves>
#include <gbs/bscanalysis.h>
namespace gbs
{
    template <typename T, size_t dim>
    class msh_edge
    {
        CurveReparametrized<T,dim> m_crv;
        BSCfunction<T> m_law;
        points_vector<T, dim> m_pts;
        T l_;
        size_t n_;
    public:
        msh_edge(const std::shared_ptr<Curve<T, dim>> &p_crv) : 
            m_crv{p_crv, abs_curv(*p_crv)},
            n_{2},
            m_law{
                std::vector<T>{ 0., 1.},
                std::vector<T>{0., 0., 1., 1.},
                1
            },
            l_{ m_crv.bounds()[1] }
        {
        }

        msh_edge(std::shared_ptr<Curve<T, dim>> &p_crv,std::array<T,2> bounds) : 
            m_crv{p_crv, abs_curv(*p_crv, bounds[0], bounds[1])},
            n_{2},
            m_law{
                std::vector<T>{ 0., 1.},
                std::vector<T>{0., 0., 1., 1.},
                1
            },
            l_{ m_crv.bounds()[1] }
        {
        }

        auto set_n_points(size_t n) -> void
        {
            n_ = fmax(2, n);
            m_pts = points_vector<T, dim>(n_);
        }

        auto update_law(T h1, T h2) -> void
        {
            points_vector<T, 1> m = {{0.}, {h1 / l_}, {1. - h2 / l_}, {1.}};
            std::vector<T> ksi = {0., 1 / (n_ - 1.), 1. - 1 / (n_ - 1.), 1.};
            m_law = interpolate(m, ksi, 2);
        }

        auto update_law(T h1, T h2, T a1, T a2) -> void
        {
            points_vector<T, 1> m = {
                {0.},
                {h1 / l_},
                {h1 / l_ + a1 * h1 / l_},
                {1. - h2 / l_ - a2 * h2 / l_},
                {1. - h2 / l_},
                {1.}};
            std::vector<T> ksi = {
                0.,
                1 / (n_ - 1.),
                2 / (n_ - 1.),
                1. - 2 / (n_ - 1.),
                1. - 1 / (n_ - 1.),
                1.};
            m_law = interpolate(m, ksi, 2);
        }

        auto compute_pnts(std::execution::parallel_policy execution = std::execution::par)
        {
            auto ksi = make_range(0., 1., n_);
            m_pts.front() = m_crv.begin();
            m_pts.back()  = m_crv.end();
            std::transform(
                // execution,
                std::next(ksi.begin()),
                std::next(ksi.end(),-1),
                std::next(m_pts.begin()),
                [&](const auto &ksi_) {
                    auto m = m_law(ksi_)* l_;
                    return m_crv(m);
                });
        }

        auto points() const noexcept -> const points_vector<T, dim> &
        {
            return m_pts;
        }
    };
}