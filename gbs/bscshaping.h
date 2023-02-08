#include <gbs/curves>
#include <gbs/bscapprox.h>

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
        std::array<Eigen::VectorXd, dim> DX;
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

        return BSCurve<T,dim>{poles,k,p};

    }

    template<typename T, size_t dim>
    auto move_to_point(BSCurve<T,dim> &crv, const point<T,dim> &pt, T u)
    {       
        crv.setPoles( moved_to_point(crv, pt, u).poles());
    }

    template<typename T, size_t dim>
    auto move_to_constraints_delta(points_vector<T, dim> &poles,const std::vector<bsc_constraint<T,dim>> &constraints_delta, const std::vector<T> &k, size_t p, size_t i1, size_t i2)
    {
        auto m = constraints_delta.size();
        auto n = i2 - i1 + 1;
        assert(m <= n );
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> B(n, m);
        for(size_t j{}; j< m; j++)
        {
            auto u = std::get<0>(constraints_delta[j]);
            auto d = std::get<2>(constraints_delta[j]);
            for (size_t i{i1}; i <= i2; i++)
            {
                B(i - i1, j) = basis_function(u, i, p, d, k);
            }
        }

        std::array<Eigen::VectorXd, dim> DX;
        std::array<Eigen::VectorXd, dim> DP;
        for(size_t i{}; i < dim; i++){
            DX[i] = Eigen::VectorXd(m);
            for(size_t j{}; j< m; j++)
            {
                auto u = std::get<0>(constraints_delta[j]);
                auto d = std::get<2>(constraints_delta[j]);
                DX[i](j) = std::get<1>(constraints_delta[j])[i];
            }
        }

        if(m < n)
        {
            auto TB = B.transpose();
            auto Inv = (TB * B).llt();
            for(size_t i{}; i < dim; i++)
            {
                DP[i] =  B * Inv.solve(DX[i]) ;
            }
        }
        else{
            auto Inv = B.partialPivLu();
            for(size_t i{}; i < dim; i++)
            {
                DP[i] =  Inv.solve(DX[i]) ;
            }
        }

        for (size_t i{i1}; i <= i2; i++)
        {
            for(size_t j{}; j < dim; j++)
            {
                poles[i][j] = poles[i][j] + DP[j](i-i1);
            }
        }
    }

    template<typename T, size_t dim>
    auto moved_to_constraints(const BSCurve<T,dim> &crv, const std::vector<bsc_constraint<T,dim>> &constraints)
    {       
        auto poles = crv.poles();
        auto n = poles.size();
        auto k     = crv.knotsFlats();
        auto p     = crv.degree();

        size_t i1{}, i2{n-1};
        std::vector<bsc_constraint<T,dim>> constraints_delta;
        for(auto &cstr : constraints)
        {
            auto u_ = std::get<0>(cstr);
            auto d_ = std::get<2>(cstr);
            auto [i1_, i2_] = find_span_range(n, p , u_ , k);
            i1 = std::min(i1,i1_);
            i2 = std::max(i2,i2_);
            constraints_delta.push_back( {u_,  std::get<1>(cstr) - crv(u_, d_), d_} );
        }
        
        move_to_constraints_delta(poles, constraints_delta, k, p , i1, i2);

        return BSCurve<T,dim>{poles,k,p};

    }

}