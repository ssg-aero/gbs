#pragma once
#include <gbs/bscurve.h>
#include <gbs/bscanalysis.h>
#include <gbs/bscinterp.h>
#include <gbs/vecop.h>
#include <gbs/transform.h>
#include <numbers>
#include <optional>
#include <type_traits>
namespace gbs

{
    /**
     * @brief check if curve's menber fullfill bspline definition
     **/ 
    template <typename T, size_t dim,bool rational>
    auto check_curve(const BSCurveGeneral<T, dim,rational> &crv)
    {
        return check_curve(crv.poles().size(),crv.knotsFlats(),crv.degree());
    }
    /**
     * @brief Set all bspline curves to the same degree
     * 
     * @tparam Container 
     * @param bs_lst 
     */
    template <typename Container>
    auto unify_curves_degree(Container &bs_lst) -> void
    {
        // Find the maximum degree among all curves
        auto max_degree = std::max_element(
            bs_lst.begin(),
            bs_lst.end(),
            [](const auto &C1, const auto &C2) {
                return C1.degree() < C2.degree();
            })->degree();

        // Increase the degree of each curve to the maximum
        std::for_each(
            bs_lst.begin(),
            bs_lst.end(),
            [max_degree](auto &C) {
                C.increaseDegree(max_degree-C.degree());
            });
    }

    // Suspicious, I guess it contains error
    // but yet required for compilation -> invetigate and delete
    template <typename T, size_t dim>
    auto unify_knots( 
              std::pair<points_vector<T, dim >, std::vector<T> > &pk1,
        const std::pair<points_vector<T, dim >, std::vector<T> > &pk2,
        size_t degree
    )
    {
        auto km2 = unflat_knots(pk2.second);
        std::for_each(
            // std::execution::par,
            km2.begin(), km2.end(),
            [&pk1, degree](const auto &km)
            {
                auto u = km.first;
                auto m = km.second - multiplicity(pk1.second, u);
                for (auto i = 0; i < m; i++)
                    insert_knot(u, degree, pk1.second, pk1.first);
            });
    }

/**
 * @brief Unifies the knots of a set of B-spline curves.
 * 
 * Given a range of B-spline curves, this function ensures that all curves 
 * have a common knot vector. The curves in the range must have the same degree.
 * 
 * @tparam It      Type of the iterator pointing to the curves.
 * @param bs_begin Iterator pointing to the beginning of the curve range.
 * @param bs_end   Iterator pointing to the end of the curve range.
 * 
 * @throw std::invalid_argument If the curves don't have the same degree.
 */
    template <typename It>
    auto unify_curves_knots(It bs_begin, It bs_end) -> void
    {
        // Check that all curves have the same degree
        auto p = bs_begin->degree();
        std::for_each(bs_begin, bs_end,[&p](const auto &C_)
        {
            if (C_.degree() != p)
                throw std::invalid_argument("unify_curves_knots: need curves with same degree");
        });

        // Compute the desired knot vector
        auto &C = *bs_begin;
        auto km = unflat_knots(C.knotsFlats());
        // Make the first curve compatible with the others
        std::for_each(
            std::next(bs_begin),
            bs_end,
            [&](auto &C_) {
                C_.changeBounds(C.bounds());
                auto km_= unflat_knots(C_.knotsFlats());
                C.insertKnots( km_ );
            }
        );

        // Apply the common knots vector to all remaining curves
        km = unflat_knots(C.knotsFlats()); 
        std::for_each(
            std::next(bs_begin),
            bs_end,
            [&](auto &C_) {
                C_.insertKnots( km );
            });
    }

    /**
     * @brief unify the knots of a bspline curve set
     * 
     * @tparam Container 
     * @param bs_lst 
     */
    template <typename Container>
    auto unify_curves_knots(Container &bs_lst) -> void
    {
        unify_curves_knots(bs_lst.begin(), bs_lst.end());
    }

    template <typename Container>
    auto unified_curves(const Container &bs_lst)
    {
        Container bs_lst_copy(bs_lst.size());
        std::transform(
            bs_lst.begin(),bs_lst.end(),
            bs_lst_copy.begin(),
            [](const auto &bs)
            {
                auto bs_copy {bs};
                return bs_copy;
            }
        );
        unify_curves_degree(bs_lst_copy);
        unify_curves_knots(bs_lst_copy);
        return bs_lst_copy;
    }
    /**
     * @brief join 2 curves at tail/head to build a new curve, geometrical definition of both curves is preserved
     * 
     * @tparam T 
     * @tparam dim 
     * @tparam rational1 
     * @tparam rational2 
     * @param crv1 
     * @param crv2 
     * @return auto 
     */
    template <typename T, size_t dim, bool rational1, bool rational2>
    auto join(const gbs::BSCurveGeneral<T, dim, rational1> &crv1,
              const gbs::BSCurveGeneral<T, dim, rational2> &crv2)
    {
        // make curves copies uniform in definition
        using crvType = typename std::conditional<rational1 || rational2,gbs::BSCurveRational<T,dim>,gbs::BSCurve<T,dim>>::type;
        //put both curves at same def
        std::vector<crvType> lst = {crvType(crv1),crvType(crv2)}; // create a cpy
        unify_curves_degree(lst);
        auto crv1_cp = lst.front();
        auto crv2_cp = lst.back(); 
        // recover data of same dimension
        auto poles1 = crv1_cp.poles();
        auto poles2 = crv2_cp.poles();
        auto k1(crv1_cp.knotsFlats());
        auto k2(crv2_cp.knotsFlats());
        auto p1 = crv1_cp.degree();
        auto p2 = crv2_cp.degree();

        // join poles
        if(rational1 || rational2)
        {
            //set tail/head with a weight of 1
            scale_poles(poles1,1. / poles1.back().back());
            scale_poles(poles2,1. / poles2.front().back());
            auto pt = T(0.5) * (weight_projection(poles1.back())+weight_projection(poles2.front()));
            poles1.back() = add_weight(pt,poles1.back().back());
        }
        else
        {
            poles1.back() = T(0.5) * (poles1.back() + poles2.front());
        }
        poles1.insert(poles1.end(), std::next(poles2.begin()), poles2.end());
        // join knots
        auto k1_end = k1.back();
        auto k2_start = k2.front();
        std::transform(
            k2.begin(),
            k2.end(),
            k2.begin(),
            [&k2_start,&k1_end](const auto &k_)
            {
                return k_ + k1_end - k2_start;
            }
        );
        // k1.pop_back();
        k1.erase( k1.end() - 1, k1.end() );
        k1.insert(k1.end() , std::next(k2.begin(),p2+1), k2.end());
        // create result
        return crvType(poles1,k1,p1);
    }

    template <typename T, size_t dim, bool rational1, bool rational2>
    auto join(const gbs::BSCurveGeneral<T, dim, rational1> *crv1,
              const gbs::BSCurveGeneral<T, dim, rational2> *crv2)// -> std::unique_ptr<gbs::BSCurveGeneral<T, dim, rational1 || rational2>>
    {
        using crvType = typename std::conditional<rational1 || rational2,gbs::BSCurveRational<T,dim>,gbs::BSCurve<T,dim>>::type;
        return std::make_unique<crvType>(join(*crv1,*crv2));
    }


    template <typename T, size_t dim>
    auto c2_connect(const Curve<T, dim> &crv1,
                    const Curve<T, dim> &crv2) -> BSCurve<T, dim>
    {
        auto u1 = crv1.bounds()[1];
        auto u2 = crv2.bounds()[0];
        std::vector<gbs::constrType<T, dim, 3>> Q =
            {
                {crv1.value(u1), crv1.value(u1, 1), crv1.value(u1, 2)},
                {crv2.value(u2), crv2.value(u2, 1), crv2.value(u2, 2)}};

        return interpolate(Q, gbs::KnotsCalcMode::CHORD_LENGTH);
    }



    template <typename T, size_t dim>
    auto c2_connect(const Curve<T, dim> &crv1,
                    const Curve<T, dim> &crv2,
                    T t1 , T t2) -> BSCurve<T, dim>
    {
        auto u1 = crv1.bounds()[1];
        auto u2 = crv2.bounds()[0];
        std::vector<gbs::constrType<T, dim, 4>> Q =
            {
                {crv1.value(u1), crv1.value(u1, 1), crv1.value(u1, 2), t1 * crv1.value(u1, 3)},
                {crv2.value(u2), crv2.value(u2, 1), crv2.value(u2, 2), t2 * crv2.value(u2, 3)}};

        return interpolate(Q, gbs::KnotsCalcMode::CHORD_LENGTH);
    }


    template <typename T>
    auto c2_connect(const point<T, 2> &p1,
                    const point<T, 2> &p2,
                    const point<T, 2> &t1,
                    const point<T, 2> &t2,
                    const point<T, 2> &c1,
                    const point<T, 2> &c2,
                    T e
                     ) -> BSCurve<T, 2>
    {
        auto d  = norm(p1 - p2);
        auto u1 = 0.;
        auto u2 = std::numbers::pi * d * 0.5;
        std::vector<bsc_constraint<T, 2>> cstr_lst = {
            bsc_constraint<T, 2>{u1,t1,1}
            ,
            bsc_constraint<T, 2>{u1,c1,2}
            ,
            bsc_constraint<T, 2>{u2,t2,1}
            ,
            bsc_constraint<T, 2>{u2,c2,2}
            ,
            // bsc_constraint<T, 2>{u2,d2,3}
            // ,
            // bsc_constraint<T, 2>{u1,d1,3}
            // ,
        };

        auto c = interpolate(
            bsc_bound<T, 2>{u1, p1}, 
            bsc_bound<T, 2>{u2, p2}, 
            cstr_lst, 
            // 5
            3
        );

        auto u_mid =0.5*(u1+u2); 
        auto cr = c(u_mid,2);
        cr = cr / norm(cr) * 0.5 * d * e;
        auto tg = c(u_mid,1);
        tg = tg / norm(tg) * 0.5 * d ;
        cstr_lst.push_back(bsc_constraint<T, 2>{u_mid,tg,1});
        cstr_lst.push_back(bsc_constraint<T, 2>{u_mid,cr,2});
        // cstr_lst.push_back(bsc_constraint<T, 2>{u1,{0.,0.},3});
        // cstr_lst.push_back(bsc_constraint<T, 2>{u2,{0.,0.},3});
        return interpolate(
            bsc_bound<T, 2>{u1, p1}, 
            bsc_bound<T, 2>{u2, p2}, 
            cstr_lst, 
            // 6
            4
        );

    }
    /**
     * @brief C2 connection mimicing ellipse but using bspline
     * 
     * @tparam T 
     * @param crv1 
     * @param crv2 
     * @param e 
     * @return BSCurve<T, 2> 
     */
    template <typename T>
    auto c2_connect(const Curve<T, 2> &crv1,
                    const Curve<T, 2> &crv2, T e) -> BSCurve<T, 2>
    {
        // TODO check if p1 == p2
        auto p1 = crv1.end();
        auto p2 = crv2.begin();
        auto t1 = crv1.end(1);
        auto t2 = crv2.begin(1);
        // auto a12= std::asin(norm(t1 ^ t2)/(norm(t1)* norm(t2)));
        auto c1 = crv1.end(2);
        auto c2 = crv2.begin(2);
        // auto d1 = crv1.end(3);
        // auto d2 = crv2.begin(3);
        return c2_connect(p1,p2,t1,t2,c1,c2,e);
    }

    template <typename T>
    auto c2_connect(const Curve<T, 2> &crv1,
                    const Curve<T, 2> &crv2, 
                    T u1, T u2,
                    bool side1, bool side2,
                    T e
                    ) -> BSCurve<T, 2>
    {
        // TODO check if p1 == p2
        T s1 = side1 ? 1. : -1.;
        T s2 = side2 ? 1. : -1.;
        auto p1 = crv1(u1);
        auto p2 = crv2(u2);
        auto t1 = s1*crv1(u1,1);
        auto t2 = s2*crv2(u2,1);
        auto c1 = s1*crv1(u1,2);
        auto c2 = s2*crv2(u2,2);

        return c2_connect(p1,p2,t1,t2,c1,c2,e);
    }

    template <typename T, size_t dim>
    auto c2_connect(const Curve<T, dim> &crv1,
                    const Curve<T, dim> &crv2,
                    T e) -> BSCurve<T, dim>
    {
        auto u1 = crv1.bounds()[1];
        auto u2 = crv2.bounds()[0];
        auto p1 = crv1.value(u1);
        auto p2 = crv2.value(u2);
        auto t1 = crv1.value(u1,1);
        auto t2 = crv2.value(u2,1);
        auto t  = (t1 - t2);
        t = t / norm(t);
        auto d = norm(p1 - p2);
        auto p = 0.5 * (p1 + p2) + e * d * t;

        auto u = gbs::curve_parametrization(std::vector<point<T,dim>>{p1,p,p2}, KnotsCalcMode::CHORD_LENGTH,true);
        std::vector<constrPoint<T, dim>> Q =
        {
            {p1, 0, 0},
            {t1, 0, 1},
            {crv1.value(u1,2), 0, 2},
            {p , 0.5*d, 0},
            {p2, d, 0},
            {t2, d, 1},
            {crv2.value(u2,2), d, 2},
        };

        T deg = 3; // degree

        std::vector<T> knotsFlats = { 0., 0. ,0.,0. , 0.25*d,  0.5*d,  0.75*d, d ,d ,d ,d};
        auto poles = gbs::build_poles(Q,knotsFlats,deg);
        return BSCurve<T,dim>(poles,knotsFlats,deg);
    }
    template <typename T, size_t dim>
    auto c3_connect(const Curve<T, dim> &crv1,
                    const Curve<T, dim> &crv2) -> BSCurve<T, dim>
    {
        auto u1 = crv1.bounds()[1];
        auto u2 = crv2.bounds()[0];
        std::vector<gbs::constrType<T, dim, 4>> Q =
            {
                {crv1.value(u1), crv1.value(u1, 1), crv1.value(u1, 2), crv1.value(u1, 3)},
                {crv2.value(u2), crv2.value(u2, 1), crv2.value(u2, 2), crv2.value(u2, 3)}};

        return interpolate(Q, gbs::KnotsCalcMode::CHORD_LENGTH);
    }

    /**
     * @brief C3 connection mimicing ellipse but using bspline
     * 
     * @tparam T 
     * @tparam dim
     * @param crv1 
     * @param crv2 
     * @param e 
     * @return BSCurve<T, 2> 
     */
    template <typename T, size_t dim>
    auto c3_connect(const Curve<T, dim> &crv1,
                    const Curve<T, dim> &crv2, T e) -> BSCurve<T, 2>
    {
        // TODO check if p1 == p2
        auto p1 = crv1.end();
        auto p2 = crv2.begin();
        auto t1 = crv1.end(1);
        auto t2 = crv2.begin(1);
        // auto a12= std::asin(norm(t1 ^ t2)/(norm(t1)* norm(t2)));
        auto c1 = crv1.end(2);
        auto c2 = crv2.begin(2);
        auto d1 = crv1.end(3);
        auto d2 = crv2.begin(3);
        auto d  = norm(p1 - p2);
        auto u1 = 0.;
        auto u2 = std::numbers::pi * d * 0.5;
        std::vector<bsc_constraint<T, 2>> cstr_lst = {
            bsc_constraint<T, 2>{u1,t1,1}
            ,
            bsc_constraint<T, 2>{u1,c1,2}
            ,
            bsc_constraint<T, 2>{u2,t2,1}
            ,
            bsc_constraint<T, 2>{u2,c2,2}
            ,
            bsc_constraint<T, 2>{u2,d2,3}
            ,
            bsc_constraint<T, 2>{u1,d1,3}
            ,
        };

        auto c = interpolate(
            bsc_bound<T, 2>{u1, p1}, 
            bsc_bound<T, 2>{u2, p2}, 
            cstr_lst, 
            5
            // 3
        );

        auto u_mid =0.5*(u1+u2); 
        auto cr = c(u_mid,2);
        cr = cr / norm(cr) * 0.5 * d * e;
        auto tg = c(u_mid,1);
        tg = tg / norm(tg) * 0.5 * d ;
        cstr_lst.push_back(bsc_constraint<T, 2>{u_mid,tg,1});
        cstr_lst.push_back(bsc_constraint<T, 2>{u_mid,cr,2});
        // cstr_lst.push_back(bsc_constraint<T, 2>{u_mid,-1.*tg,3});
        return interpolate(
            bsc_bound<T, 2>{u1, p1}, 
            bsc_bound<T, 2>{u2, p2}, 
            cstr_lst, 
            // 6
            4
        );

    }

    /**
     * @brief Create a curves's copy with and additional dimension. The value of the dimension can be specified, the default is 0.
     * 
     * @tparam T 
     * @tparam dim 
     * @tparam rational 
     * @param crv 
     * @param val 
     * @return auto 
    **/
    template <typename T, size_t dim, bool rational>
    auto add_dimension(const BSCurveGeneral<T, dim, rational> &crv,T val=0.)
    {
        points_vector<T,dim+rational+1> poles(crv.poles().size());
        std::transform(
            std::execution::par,
            crv.poles().begin(),
            crv.poles().end(),
            poles.begin(),
            [&val](const auto &p_){
                return add_dimension(p_,val);
            }
        );
        using bs_type = typename std::conditional<rational,BSCurveRational<T, dim+1>,BSCurve<T, dim+1>>::type;
        return bs_type( poles,  crv.knotsFlats(), crv.degree() );
    }
/**
 * @brief Build a curve as an extention of the curve crv to the point pt from position u folowing curve's direction
 * 
 * @tparam T 
 * @tparam dim 
 * @tparam rational 
 * @param crv    : the curve to be extended
 * @param pt     : point toward the curve has to be extended
 * @param u      : position of the begining of extention
 * @param u_new  : parameter corresponding of pt on extention
 * @param natural_end : toggle null curvature at pt
 * @param max_cont : max continuity, if nothing provided crv's degree is used
 * @return auto BSCurve<T,dim>
 */
    template <typename T, size_t dim, bool rational>
    auto extention_to_point(const BSCurveGeneral<T,dim,rational> &crv, const point<T,dim> &pt, T u, T u_new, bool natural_end , std::optional<size_t> max_cont = std::nullopt)
    { 
        using bsc_type = typename std::conditional<rational,BSCurveRational<T, dim>,BSCurve<T, dim>>::type;
        gbs::bsc_bound<T,dim> pt_begin = {u,crv(u)};
        gbs::bsc_bound<T,dim> pt_end   = {u_new,pt};

        auto p = crv.degree();
        std::vector<gbs::bsc_constraint<T,dim>> cstr_lst;
        for(int i = 1 ; (i < p) && (i <= max_cont.value_or(p)) ; i++)
        {
            cstr_lst.push_back({u,crv(u,i),i});
        }
        if(p>2 && natural_end) 
        {
            point<T,dim> cu0;
            cu0.fill(0.);
            cstr_lst.push_back({u_new,cu0,2}); //Natural BS
        }
        return bsc_type{gbs::interpolate<T,dim>(pt_begin,pt_end,cstr_lst,std::min(p,cstr_lst.size()+1))};
    }
/**
 * @brief Build a curve as an extention of the curve crv to the point pt from position u following curve's direction
 * 
 * @tparam T 
 * @tparam dim 
 * @tparam rational 
 * @param crv    : the curve to be extended
 * @param pt     : point toward the curve has to be extended
 * @param u      : position of the begining of extention
 * @param natural_end : toggle null curvature at pt
 * @param max_cont : max continuity, if nothing provided crv's degree is used
 * @return auto BSCurve<T,dim>
 */
    template <typename T, size_t dim, bool rational>
    auto extention_to_point(const BSCurveGeneral<T,dim,rational> &crv, const point<T,dim> &pt, T u, bool natural_end, std::optional<size_t> max_cont = std::nullopt)
    {
        auto t = crv(u,1);
        auto dl = gbs::norm(crv.end() - pt);
        auto du = dl / gbs::norm(t);
        auto u_new = u + dl * du;
        return extention_to_point(crv,pt,u,u_new,natural_end,max_cont);
    }
    /**
     * @brief Extend a curve to the point pt from position u following curve's direction
     * 
     * @tparam T 
     * @tparam dim 
     * @tparam rational 
     * @param crv    : the curve to be extended
     * @param pt     : point toward the curve has to be extended
     * @param natural_end : toggle null curvature at pt
     * @param max_cont : max continuity, if nothing provided crv's degree is used
     * @return auto BSCurve<T,dim>
     */
    template <typename T, size_t dim, bool rational>
    auto extended_to_point(const BSCurveGeneral<T,dim,rational> &crv, const point<T,dim> &pt, bool natural_end, std::optional<size_t> max_cont = std::nullopt)
    {
        auto u = crv.bounds()[1];
        auto extention =extention_to_point(crv,pt,u,natural_end,max_cont); 
        return join( crv, extention );
    }
    /**
     * @brief Extend curve's end/start to point, the curve definition is kept, i.e. between original curve's bound definition are identical
     * 
     * @tparam T 
     * @tparam dim 
     * @param crv 
     * @param pt 
     * @param at_end 
     * @param natural_end : toggle null curvature at pt
     * @param max_cont : max continuity, if nothing provided crv's degree is used
     * @return BSCurve<T,dim> 
     */
    template <typename T, size_t dim, bool rational>
    auto extended_to_point(const BSCurveGeneral<T,dim,rational> &crv, const point<T,dim> &pt, bool at_end, bool natural_end, std::optional<size_t> max_cont = std::nullopt)
    {
        if(at_end) return extended_to_point(crv,pt,natural_end,max_cont);
        using bs_type = typename std::conditional<rational,BSCurveRational<T, dim>,BSCurve<T, dim>>::type ;
        auto crv_rev = bs_type{crv};
        auto [u1, u2] = crv_rev.bounds();
        crv_rev.reverse();
        crv_rev = extended_to_point(crv_rev,pt,natural_end,max_cont);
        crv_rev.reverse();
        // reset parametrization as the original curve
        auto du = crv_rev.bounds()[1] - u2;
        auto k = crv_rev.knotsFlats();
        std::transform(k.begin(),k.end(),k.begin(),[du](const auto k_){return k_-du;});
        return bs_type{crv_rev.poles(),k,crv.degree()};
    }
    /**
     * @brief Extend curve's end by length, the curve definition is kept, i.e. between original curve's bound definition are identical
     * 
     * @tparam T 
     * @tparam dim 
     * @param crv 
     * @param l 
     * @param at_end 
     * @param relative 
     * @param natural_end : toggle null curvature at pt
     * @param max_cont : max continuity, if nothing provided crv's degree is used
     * @return auto 
     */
    template <typename T, size_t dim, bool rational>
    auto extended(const BSCurveGeneral<T,dim,rational> &crv, T l, bool at_end, bool relative, bool natural_end, std::optional<size_t> max_cont = std::nullopt)
    {
        point<T,dim> t,p;
        if(relative)
        {
            l *= length(crv);
        }
        if(at_end)
        {
            t  = crv.end(1);
            t = t / norm(t);
            p = crv.end();
        }
        else
        {
            t  = -1.*crv.begin(1);
            t = t / norm(t);
            p = crv.begin();
        }
        
        return extended_to_point(crv,p+t*l,at_end,natural_end,max_cont);
        
    }

} // namespace gbs