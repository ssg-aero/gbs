#include <gbs/curves>

namespace gbs
{
    template <typename T, size_t dim>
    auto move_to_point(points_vector<T, dim> &poles, const point<T, dim> &D, const std::vector<T> &k, size_t p, T u, size_t i1, size_t i2)
    {
        auto n = i2 - i1 + 1;
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> B(n, 1);
        for (size_t i{i1}; i <= i2; i++)
        {
            B(i - i1, 0) = basis_function(u, i, p, 0, k);
        }
        auto TB = B.transpose();
        auto Inv = (TB * B).llt();
        Eigen::VectorXd Dx(1), Dy(1), Dz(1);
        Dx(0) = D[0];
        Dy(0) = D[1];
        Dz(0) = D[2];

        auto dPx = Inv.solve(Dx) * TB;
        auto dPy = Inv.solve(Dy) * TB;
        auto dPz = Inv.solve(Dz) * TB;

        for (size_t i{i1}; i <= i2; i++)
        {
            poles[i] = poles[i] + point<T, 3>{dPx(i - i1), dPy(i - i1), dPz(i - i1)};
        }

    }

    template<typename T, size_t dim>
    auto moved_to_point(const BSCurve<T,dim> &crv, const point<T,dim> &pt, T u)
    {       
        auto poles = crv.poles();
        auto n = poles.size();
        auto k     = crv.knotsFlats();
        auto p     = crv.degree();

        auto [i1, i2] = find_span_range(n, p , u , k);
        auto D = pt - crv(u);

        move_to_point(poles, D, k, p , u, i1, i2);

        return BSCurve<T,3>{poles,k,p};

    }

}