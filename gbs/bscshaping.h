#include <gbs/curves>
#include <gbs/bscapprox.h>

namespace gbs
{

    /**
      @brief Edit poles to pass through constraints
     * 
     * @tparam T 
     * @tparam dim 
     * @param poles 
     * @param constraints_delta 
     * @param k 
     * @param p 
     * @param i1 
     * @param i2 
     * @return auto 
     */
    template<typename T, size_t dim>
    auto move_to_constraints_delta(points_vector<T, dim> &poles,const std::vector<bsc_constraint<T,dim>> &constraints_delta, const std::vector<T> &k, size_t p, size_t i1, size_t i2)
    {
        auto m = constraints_delta.size();
        auto n = i2 - i1 + 1;
 
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> B(m, n);
        for(size_t i{}; i< m; i++)
        {
            auto u = std::get<0>(constraints_delta[i]);
            auto d = std::get<2>(constraints_delta[i]);
            for (size_t j{i1}; j <= i2; j++)
            {
                B(i, j - i1) = basis_function(u, j, p, d, k);
            }
        }

        std::array<Eigen::VectorXd, dim> DX;
        std::array<Eigen::VectorXd, dim> DP;
        for(size_t i{}; i < dim; i++){
            DX[i] = Eigen::VectorXd(m);
            for(size_t j{}; j< m; j++)
            {
                DX[i](j) = std::get<1>(constraints_delta[j])[i];
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

        for (size_t i{i1}; i <= i2; i++)
        {
            for(size_t j{}; j < dim; j++)
            {
                poles[i][j] = poles[i][j] + DP[j](i-i1);
            }
        }
    }


    template<typename T, size_t dim>
    auto moved_to_constraints_delta(const BSCurve<T,dim> &crv, const std::vector<bsc_constraint<T,dim>> &constraints_delta, size_t i1, size_t i2)
    {
        auto poles = crv.poles();
        auto n = poles.size();
        auto k     = crv.knotsFlats();
        auto p     = crv.degree();
        move_to_constraints_delta(poles, constraints_delta, k, p , i1, i2);
        return BSCurve<T,dim>{poles,k,p};
    }

    template<typename T, size_t dim> 
    auto find_span_range(const BSCurve<T,dim> &crv, const std::vector<bsc_constraint<T,dim>> &constraints)
    {
        auto poles = crv.poles();
        auto n = poles.size();
        auto k     = crv.knotsFlats();
        auto p     = crv.degree();
        size_t i1{n-1}, i2{};
        for(auto &cstr : constraints)
        {
            auto u_ = std::get<0>(cstr);
            auto [i1_, i2_] = find_span_range(n, p , u_ , k);
            i1 = std::min(i1,i1_);
            i2 = std::max(i2,i2_);
        }
        return std::make_tuple(i1, i2);
    }

    template<typename T, size_t dim>
    auto convert_constraints_to_delta(const BSCurve<T,dim> &crv, const std::vector<bsc_constraint<T,dim>> &constraints)
    {
        std::vector<bsc_constraint<T,dim>> constraints_delta;
        for(auto &cstr : constraints)
        {
            auto u_ = std::get<0>(cstr);
            auto d_ = std::get<2>(cstr);
            constraints_delta.push_back( {u_,  std::get<1>(cstr) - crv(u_, d_), d_} );
        }
        return constraints_delta;
    }

    /**
     * @brief Builds crv copy and edit poles to match constraints
     * 
     * @tparam T 
     * @tparam dim 
     * @param crv 
     * @param constraints 
     * @return auto 
     */
    template<typename T, size_t dim>
    auto moved_to_constraints(const BSCurve<T,dim> &crv, const std::vector<bsc_constraint<T,dim>> &constraints, bool fix_bounds=false)
    {       
        auto constraints_delta = convert_constraints_to_delta(crv, constraints);
        auto n = crv.poles().size();
        auto [i1, i2] = find_span_range(crv, constraints);
        if(fix_bounds)
        {
            i1 = std::max<size_t>(1,i1);
            i2 = std::min<size_t>(n-2,i2);
        }

        return moved_to_constraints_delta(crv, constraints, i1, i2);
    }

    template<typename T, size_t dim>
    auto move_to_constraints(BSCurve<T,dim> &crv, const std::vector<bsc_constraint<T,dim>> &constraints, bool fix_bounds=false)
    {       
        crv.copyPoles( moved_to_constraints(crv, constraints, fix_bounds).poles() );
    }
    /**
     * @brief Edit poles to pass through pt at parameter u
     * 
     * @tparam T 
     * @tparam dim 
     * @param poles 
     * @param D 
     * @param k 
     * @param p 
     * @param u 
     * @param i1 
     * @param i2 
     * @return auto 
     */
    template <typename T, size_t dim>
    auto move_to_point_delta(points_vector<T, dim> &poles, const point<T, dim> &D, const std::vector<T> &k, size_t p, T u, size_t i1, size_t i2)
    {
        move_to_constraints_delta(poles, {{u,D,0}}, k, p, i1, i2);
    }

    /**
     * @brief Builds crv copy and edit poles to pass through pt at parameter u
     * 
     * @tparam T 
     * @tparam dim 
     * @param crv 
     * @param pt 
     * @param u 
     * @param fix_bounds
     * @return BSCurve<T,dim> 
     */
    template<typename T, size_t dim>
    auto moved_to_point(const BSCurve<T,dim> &crv, const point<T,dim> &pt, T u, bool fix_bounds=false)
    {       
        auto poles = crv.poles();
        auto n = poles.size();
        auto k     = crv.knotsFlats();
        auto p     = crv.degree();

        auto [i1, i2] = find_span_range(n, p , u , k);
        if(fix_bounds)
        {
            i1 = std::max<size_t>(1,i1);
            i2 = std::min<size_t>(n-2,i2);
        }
        auto D = pt - crv(u);

        move_to_point_delta(poles, D, k, p , u, i1, i2);

        return BSCurve<T,dim>{poles,k,p};

    }

    /**
     * @brief Edit curve's poles to pass through pt at parameter u
     * 
     * @tparam T 
     * @tparam dim 
     * @param crv 
     * @param pt 
     * @param u 
     * @param fix_bounds
     * @return void 
     */
    template<typename T, size_t dim>
    auto move_to_point(BSCurve<T,dim> &crv, const point<T,dim> &pt, T u, bool fix_bounds=false)
    {       
        crv.copyPoles( moved_to_point(crv, pt, u, fix_bounds).poles() );
    }

}