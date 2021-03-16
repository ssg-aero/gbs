#include <gbs/bscurve.h>
namespace gbs
{
    template <typename T, size_t dim>
    class msh_edge
    {
        std::shared_ptr<Curve<T, dim>> p_crv_;
        std::array<T, 2> bounds_;
        BSCurve<T, 1> law_;
        BSCurve<T, 1> abs_crv_law_;
        points_vector<T, dim> pts_;
        T l_;
        size_t n_;

    public:
        msh_edge(Curve<T, dim> *p_crv, std::array<T, 2> bounds = {0., 1.}) : p_crv_{p_crv}, bounds_{bounds}, n_{2}
        {
            abs_crv_law_ = abs_curv(*p_crv);
            law_ = BSCurve<T, 1>{
                points_vector<T, 1>{
                    {{bounds[0]},
                     {bounds[1]}}},
                std::vector<T>{
                    0., 0., 1., 1.},
                1};
            l_ = length(*p_crv);
        }

        auto set_points(size_t n) -> void
        {
            n_ = fmax(2, n);
            pts_ = points_vector<T, dim>(n);
        }

        auto update_law(T h1, T h2) -> void
        {
            points_vector<T, 1> m = {{0.}, {h1 / l_}, {1. - h2 / l_}, {1.}};
            std::vector<T> ksi = {0., 1 / (n_ - 1.), 1. - 1 / (n_ - 1.), 1.};
            law_ = interpolate(m, ksi, 2);
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
            law_ = interpolate(m, ksi, 2);
        }

        auto compute_pnts(std::execution::parallel_policy execution = std::execution::par)
        {
            auto ksi = make_range(0., 1., n_);
            std::transform(
                execution,
                ksi.begin(),
                ksi.end(),
                pts_.begin(),
                [&](const auto &ksi_) {
                    auto m = law_.value(ksi_)[0] * l_;
                    auto u = abs_crv_law_(m)[0];
                    return p_crv_->value(u);
                });
        }

        auto points() const noexcept -> const points_vector<T, dim> &
        {
            return pts_;
        }
    };
}