#pragma once
#include <gbs/bscurve.h>
#include <gbs/bscinterp.h>
#include <gbs/vecop.h>
#include <gbs/transform.h>
#include <numbers>

namespace gbs

{
    /**
     * @brief check if curve's menber fullfill bspline definition
     **/ 
    template <typename T, size_t dim,bool rational>
    auto check_curve(const BSCurveGeneral<T, dim,rational> &crv)
    {
        return check_curve(crv.poles(),crv.knotsFlats(),crv.degree());
    }
    /**
     * @brief Set all bspline curves to the same degree
     * 
     * @tparam Container 
     * @param bs_lst 
     */
    template <typename Container>
    auto unify_degree(Container &bs_lst) -> void
    {
        auto &C = bs_lst.front();
        // Set first curve at maximun degree
        std::for_each(
            std::next(bs_lst.begin()),
            bs_lst.end(),
            [&](auto &C_) {
                while (C.degree()< C_.degree())
                {
                    C.increaseDegree();
                } 
            });
        // Set the others at first curve's degree
        std::for_each(
            std::next(bs_lst.begin()),
            bs_lst.end(),
            [&](auto &C_) {
                while (C_.degree()< C.degree())
                {
                    C_.increaseDegree();
                } 
            });
    }
    /**
     * @brief unify the knots of a bspline curve set
     * 
     * @tparam Container 
     * @param bs_lst 
     */
    template <typename Container>
    auto unify_knots(Container &bs_lst) -> void
    {
        auto p = bs_lst.front().degree();
        std::for_each(bs_lst.begin(),bs_lst.end(),[&p](const auto &C_)
        {
            if (C_.degree() != p)
                throw std::exception("unify_knots: need curves with same degree");
        });

        auto &C = bs_lst.front();

        std::for_each(
            std::next(bs_lst.begin()),
            bs_lst.end(),
            [&](auto &C_) {
                C_.changeBounds(C.bounds());
                auto km_ = unflat_knots(C_.knotsFlats());
                C.insertKnots( km_ );
            });

        auto km = unflat_knots(C.knotsFlats());

        std::for_each(
            std::next(bs_lst.begin()),
            bs_lst.end(),
            [&](auto &C_) {
                C_.insertKnots( km );
            });
    }

    template <typename T, size_t dim, bool rational1, bool rational2>
    auto join(const gbs::BSCurveGeneral<T, dim, rational1> *crv1,
              const gbs::BSCurveGeneral<T, dim, rational2> *crv2)// -> std::unique_ptr<gbs::BSCurveGeneral<T, dim, rational1 || rational2>>
    {
        // make curves copies uniform in definition
        typedef std::conditional<rational1 || rational2,gbs::BSCurveRational<T,dim>,gbs::BSCurve<T,dim>>::type crvType;
        std::unique_ptr< gbs::BSCurveGeneral<T, dim, rational1 || rational2> > crv1_cp,crv2_cp;
        crv1_cp = std::make_unique<crvType>(*crv1); 
        crv2_cp = std::make_unique<crvType>(*crv2);
        // recover data of same dimension
        auto poles1 = crv1_cp->poles();
        auto poles2 = crv2_cp->poles();
        auto k1(crv1_cp->knotsFlats());
        auto k2(crv2_cp->knotsFlats());
        auto p1 = crv1_cp->degree();
        auto p2 = crv2_cp->degree();      
        while (p1<p2)
        {
            increase_degree(k1,poles1,p1);
            p1++;
        }
        while (p2<p1)
        {
            increase_degree(k2,poles2,p2);
            p2++;
        }
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
        k1.erase( k1.end() - p1, k1.end() );
        k1.insert(k1.end() , std::next(k2.begin(),p2), k2.end());
        // create result
        return std::make_unique<crvType>(poles1,k1,p1);
    }

    template <typename T, size_t dim, bool rational1, bool rational2>
    auto join(const gbs::BSCurveGeneral<T, dim, rational1> &crv1,
              const gbs::BSCurveGeneral<T, dim, rational2> &crv2)// -> std::unique_ptr<gbs::BSCurveGeneral<T, dim, rational1 || rational2>>
    {
        // make curves copies uniform in definition
        typedef std::conditional<rational1 || rational2,gbs::BSCurveRational<T,dim>,gbs::BSCurve<T,dim>>::type crvType;
        //put both curves at same def
        std::vector<crvType> lst = {crvType(crv1),crvType(crv2)}; // create a cpy
        unify_degree(lst);
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
    auto c2_connect(const BSCurveGeneral<T, dim, rational1> &crv1,
                    const BSCurveGeneral<T, dim, rational2> &crv2) -> BSCurve<T, dim>
    {
        auto u1 = crv1.bounds()[1];
        auto u2 = crv2.bounds()[0];
        std::vector<gbs::constrType<T, dim, 3>> Q =
            {
                {crv1.value(u1), crv1.value(u1, 1), crv1.value(u1, 2)},
                {crv2.value(u2), crv2.value(u2, 1), crv2.value(u2, 2)}};

        return interpolate(Q, gbs::KnotsCalcMode::CHORD_LENGTH);
    }



    template <typename T, size_t dim, bool rational1, bool rational2>
    auto c2_connect(const BSCurveGeneral<T, dim, rational1> &crv1,
                    const BSCurveGeneral<T, dim, rational2> &crv2,
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
    /**
     * @brief C2 connection mimicing ellipse but using bspline
     * 
     * @tparam T 
     * @tparam rational1 
     * @tparam rational2 
     * @param crv1 
     * @param crv2 
     * @param e 
     * @return BSCurve<T, 2> 
     */
    template <typename T, bool rational1, bool rational2>
    auto c2_connect(const BSCurveGeneral<T, 2, rational1> &crv1,
                    const BSCurveGeneral<T, 2, rational2> &crv2, T e) -> BSCurve<T, 2>
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
        std::vector<bsc_constrain<T, 2>> cstr_lst = {
            bsc_constrain<T, 2>{u1,t1,1}
            ,
            bsc_constrain<T, 2>{u1,c1,2}
            ,
            bsc_constrain<T, 2>{u2,t2,1}
            ,
            bsc_constrain<T, 2>{u2,c2,2}
            ,
            // bsc_constrain<T, 2>{u2,d2,3}
            // ,
            // bsc_constrain<T, 2>{u1,d1,3}
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
        cstr_lst.push_back(bsc_constrain<T, 2>{u_mid,tg,1});
        cstr_lst.push_back(bsc_constrain<T, 2>{u_mid,cr,2});
        // cstr_lst.push_back(bsc_constrain<T, 2>{u1,{0.,0.},3});
        // cstr_lst.push_back(bsc_constrain<T, 2>{u2,{0.,0.},3});
        return interpolate(
            bsc_bound<T, 2>{u1, p1}, 
            bsc_bound<T, 2>{u2, p2}, 
            cstr_lst, 
            // 6
            4
        );

    }

    template <typename T, size_t dim, bool rational1, bool rational2>
    auto c2_connect(const BSCurveGeneral<T, dim, rational1> &crv1,
                    const BSCurveGeneral<T, dim, rational2> &crv2,
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
    template <typename T, size_t dim, bool rational1, bool rational2>
    auto c3_connect(const BSCurveGeneral<T, dim, rational1> &crv1,
                    const BSCurveGeneral<T, dim, rational2> &crv2) -> BSCurve<T, dim>
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
     * @tparam rational1 
     * @tparam rational2 
     * @param crv1 
     * @param crv2 
     * @param e 
     * @return BSCurve<T, 2> 
     */
    template <typename T, bool rational1, bool rational2>
    auto c3_connect(const BSCurveGeneral<T, 2, rational1> &crv1,
                    const BSCurveGeneral<T, 2, rational2> &crv2, T e) -> BSCurve<T, 2>
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
        std::vector<bsc_constrain<T, 2>> cstr_lst = {
            bsc_constrain<T, 2>{u1,t1,1}
            ,
            bsc_constrain<T, 2>{u1,c1,2}
            ,
            bsc_constrain<T, 2>{u2,t2,1}
            ,
            bsc_constrain<T, 2>{u2,c2,2}
            ,
            bsc_constrain<T, 2>{u2,d2,3}
            ,
            bsc_constrain<T, 2>{u1,d1,3}
            ,
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
        cstr_lst.push_back(bsc_constrain<T, 2>{u_mid,tg,1});
        cstr_lst.push_back(bsc_constrain<T, 2>{u_mid,cr,2});
        // cstr_lst.push_back(bsc_constrain<T, 2>{u_mid,-1.*tg,3});
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
        points_vector<T,dim+1> poles(crv.poles().size());
        std::transform(
            std::execution::par,
            crv.poles().begin(),
            crv.poles().end(),
            poles.begin(),
            [&val](const auto &p_){
                return add_dimension(p_,val);
            }
        );
        typedef std::conditional<rational,BSCurveRational<T, dim+1>,BSCurve<T, dim+1>>::type bs_type;
        return bs_type( poles,  crv.knotsFlats(), crv.degree() );
    }

} // namespace gbs