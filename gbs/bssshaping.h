#include <gbs/surfaces>
#include <gbs/bssapprox.h>

namespace gbs
{


    /**
     * @brief Edit poles to pass through constraints
     * 
     * @tparam T 
     * @tparam dim 
     * @param poles 
     * @param constraints_delta 
     * @param ku 
     * @param kv 
     * @param p 
     * @param q 
     * @param i1 
     * @param i2 
     * @param j1 
     * @param j2 
     * @return auto 
     */
    template<typename T, size_t dim>
    auto move_to_constraints_delta(points_vector<T, dim> &poles,const std::vector<bss_constraint<T,dim>> &constraints_delta, const std::vector<T> &ku, const std::vector<T> &kv, size_t p, size_t q, size_t i1, size_t i2, size_t j1, size_t j2)
    {
        auto m = constraints_delta.size();
        auto nu = i2 - i1 + 1;
        auto nv = j2 - j1 + 1;
        auto n = nu*nv;
        auto npoles_u = ku.size() - p - 1;
 
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> B(m, n);
        T Nu{},Nv{};
        for(size_t r{}; r< m; r++)
        {
            auto u = std::get<0>(constraints_delta[r]);
            auto v = std::get<1>(constraints_delta[r]);
            auto du= std::get<3>(constraints_delta[r]);
            auto dv= std::get<4>(constraints_delta[r]);

            for (size_t j{j1}; j <= j2; j++)
            {
                Nv = basis_function(v, j, q, 0, kv);
                for (size_t i{i1}; i <= i2; i++)
                {
                    Nu = basis_function(u, i, p, 0, ku);
                    B(r, i - i1 + nu * ( j - j1) ) = Nu * Nv;
                }
            }
        }

        std::array<Eigen::VectorXd, dim> DX;
        std::array<Eigen::VectorXd, dim> DP;
        for(size_t i{}; i < dim; i++){
            DX[i] = Eigen::VectorXd(m);
            for(size_t j{}; j< m; j++)
            {
                DX[i](j) = std::get<2>(constraints_delta[j])[i];
            }
        }

        if(m < n)
        {
            auto TB = B.transpose();
            auto Inv = (B * TB).llt();
            for(size_t i{}; i < dim; i++)
            {
                DP[i] =  TB * Inv.solve(DX[i]) ;
            }
        }
        else if(m > n)
        {
            auto TB = B.transpose();
            auto Inv = (TB * B).llt();
            for(size_t i{}; i < dim; i++)
            {
                DP[i] =   Inv.solve(TB *DX[i]) ;
            }
        }
        else{
            auto Inv = B.partialPivLu();
            for(size_t i{}; i < dim; i++)
            {
                DP[i] =  Inv.solve(DX[i]) ;
            }
        }

        for (size_t j{j1}; j <= j2; j++)
        {
            for (size_t i{i1}; i <= i2; i++)
            {
                for (size_t d{}; d < dim; d++)
                {
                    auto ip = i + npoles_u * j;
                    auto ie = i - i1 + nu * ( j - j1);
                    poles[ip][d] = poles[ip][d] + DP[d](ie);
                }
            }
        }
    }

    template<typename T, size_t dim>
    auto moved_to_constraints_delta(const BSSurface<T,dim> &srf,const std::vector<bss_constraint<T,dim>> &constraints_delta, size_t i1, size_t i2, size_t j1, size_t j2)
    {
        auto poles = srf.poles();
        auto nu = srf.nPolesU();
        auto nv = srf.nPolesV();
        const auto &ku = srf.knotsFlatsU();
        const auto &kv = srf.knotsFlatsV();
        auto p = srf.degreeU();
        auto q = srf.degreeV();

        if(constraints_delta.size()>0)
            move_to_constraints_delta(poles, constraints_delta, ku, kv, p, q, i1, i2, j1, j2);

        return BSSurface<T,dim>{poles,ku, kv, p, q};
    }

    template<typename T, size_t dim> 
    auto find_span_range(const BSSurface<T,dim> &srf, const std::vector<bss_constraint<T,dim>> &constraints)
    {
        auto poles = srf.poles();
        auto nu = srf.nPolesU();
        auto nv = srf.nPolesV();
        const auto &ku = srf.knotsFlatsU();
        const auto &kv = srf.knotsFlatsV();
        auto p = srf.degreeU();
        auto q = srf.degreeV();

        size_t i1{nu-1}, i2{};
        size_t j1{nv-1}, j2{};
        for(auto &cstr : constraints)
        {
            auto u_ = std::get<0>(cstr);
            auto v_ = std::get<1>(cstr);
            auto [i1_, i2_] = find_span_range(nu, p, u_, ku);
            auto [j1_, j2_] = find_span_range(nv, q, v_, kv);
            i1 = std::min(i1,i1_);
            i2 = std::max(i2,i2_);
            j1 = std::min(j1,j1_);
            j2 = std::max(j2,j2_);
        }
        return std::make_tuple(i1, i2, j1, j2);
    }

    template<typename T, size_t dim>
    auto convert_constraints_to_delta(const BSSurface<T,dim> &srf, const std::vector<bss_constraint<T,dim>> &constraints)
    {
        std::vector<bss_constraint<T,dim>> constraints_delta;
        for(auto &cstr : constraints)
        {
            auto u_ = std::get<0>(cstr);
            auto v_ = std::get<1>(cstr);
            auto du_= std::get<3>(cstr);
            auto dv_= std::get<4>(cstr);
            constraints_delta.push_back( {u_, v_,  std::get<2>(cstr) - srf(u_, v_, du_, dv_), du_, dv_} );
        }
        return constraints_delta;
    }

    /**
     * @brief Builds srf copy and edit poles to match constraints
     * 
     * @tparam T 
     * @tparam dim 
     * @param srf 
     * @param constraints 
     * @return auto 
     */
    template<typename T, size_t dim>
    auto moved_to_constraints(const BSSurface<T,dim> &srf, const std::vector<bss_constraint<T,dim>> &constraints, bool fix_bounds=false)
    {       

        auto nu = srf.nPolesU();
        auto nv = srf.nPolesV();

        auto [i1, i2, j1, j2] = find_span_range(srf, constraints);
        if(fix_bounds)
        {
            i1 = std::max<size_t>(1,i1);
            i2 = std::min<size_t>(nu-2,i2);
            j1 = std::max<size_t>(1,j1);
            j2 = std::min<size_t>(nv-2,j2);
        }

        auto constraints_delta = convert_constraints_to_delta(srf, constraints);

        return moved_to_constraints_delta(srf, constraints_delta, i1, i2, j1, j2);
    }
    
    /**
     * @brief 
     * 
     * @tparam T 
     * @tparam dim 
     * @param srf 
     * @param pt 
     * @param u 
     * @param v 
     * @return auto 
     */
    template <typename T, size_t dim>
    auto moved_to_point(const BSSurface<T, dim> &srf, const point<T, dim> &pt, T u, T v, bool fix_bounds=false)
    {
        return moved_to_constraints(
            srf, 
            {{u,v, pt, 0, 0}},
            fix_bounds
        );
    }

    template <typename T, size_t dim>
    auto move_to_point(BSSurface<T, dim> &srf, const point<T, dim> &pt, T u, T v, bool fix_bounds=false)
    {
        return srf.copyPoles( 
            moved_to_constraints(
                srf, 
                {{u,v, pt, 0, 0}},
                fix_bounds
            ).poles()
        );
    }
}
