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

    ///
    ///@brief Search in container if at least one curve pointer has a rational definition
    ///
    ///@tparam Ptr_Container 
    ///@param p_crv_lst 
    ///@return auto 
    ///
    template <typename T, size_t  dim>
    auto has_nurbs(const std::list<Curve<T, dim>*> &p_crv_lst)
    {
        return 
        std::find_if(
            p_crv_lst.begin(),
            p_crv_lst.end(),
            [](const auto p_c){
                auto p_nurbs = dynamic_cast<BSCurveRational<T,dim>*>(p_c);
                return  p_nurbs != nullptr;
                }
        )
        != p_crv_lst.end();
    }
    /**
     * Converts a given curve to a B-spline curve (BSCurve).
     * 
     * This function takes a shared pointer to a Curve and attempts to convert it
     * into a B-spline curve representation. If the provided curve is already a 
     * BSCurve, it is returned as-is. Otherwise, it is approximated to a BSCurve
     * using the provided deviation, degree of approximation, and mode for knot calculation.
     * 
     * @param crv A shared pointer to the Curve to be converted.
     * @param dev The maximum deviation for the curve approximation (default: 0.01).
     * @param np The number of points to be used in the approximation (default: 100).
     * @param deg_approx The degree of the B-spline curve to be used in the approximation (default: 5).
     * @param mode The method to be used for knot calculation (default: CHORD_LENGTH).
     * @return A shared pointer to the resulting BSCurve.
     */
    template <typename T, size_t dim>
    auto to_bs_curve(const std::shared_ptr < Curve<T, dim> > &crv, T dev = 0.01, size_t np = 100, size_t deg_approx = 5, gbs::KnotsCalcMode mode = gbs::KnotsCalcMode::CHORD_LENGTH) -> std::shared_ptr < BSCurve<T, dim> >
    {
        auto p_bs = std::dynamic_pointer_cast<const BSCurve<T, dim> >(crv);
        if (p_bs)
        {
            return std::make_shared<BSCurve<T, dim>>( BSCurve<T, dim>{*p_bs} );
        }
        else
        {
            return std::make_shared<BSCurve<T, dim>>( approx<T, dim>(*crv, dev, deg_approx, mode, np) );
        }
    }
    /**
     * @brief check if curve's menber fullfill bspline definition
     **/ 
    template <typename T, size_t dim,bool rational>
    auto check_curve(const BSCurveGeneral<T, dim,rational> &crv)
    {
        return check_curve(crv.poles().size(),crv.knotsFlats(),crv.degree());
    }
    /**
     * Extracts essential information from a range of B-spline curves or shared pointers to B-spline curves.
     *
     * This function iterates over a range of B-spline curves or shared pointers to B-spline curves, either rational or non-rational,
     * and extracts information from each curve into a `BSCurveInfo` object. The function handles both direct objects and shared pointers,
     * making it versatile for different use cases.
     *
     * @tparam T Floating point type used for the curve coordinates and knots.
     * @tparam dim Dimensionality of the control points (poles).
     * @tparam InputIt Iterator type, which must point to either `BSCurve<T, dim>`, `BSCurveRational<T, dim>`,
     *                 `std::shared_ptr<BSCurve<T, dim>>`, or `std::shared_ptr<BSCurveRational<T, dim>>`.
     * @param first The beginning of the range of B-spline curve elements or shared pointers to them.
     * @param last The end of the range of B-spline curve elements or shared pointers to them.
     * @return A vector of `BSCurveInfo<T, dim + (rational ? 1 : 0)>` containing the extracted information for each curve in the range.
     */
    template <std::floating_point T, size_t dim, typename InputIt>
    auto get_bs_curves_info(InputIt first, InputIt last)
    {
        static_assert(
            std::is_same_v<typename std::iterator_traits<InputIt>::value_type, BSCurve<T, dim>> ||
            std::is_same_v<typename std::iterator_traits<InputIt>::value_type, BSCurveRational<T, dim>> ||
            std::is_same_v<typename std::iterator_traits<InputIt>::value_type, std::shared_ptr<BSCurve<T, dim>>> ||
            std::is_same_v<typename std::iterator_traits<InputIt>::value_type, std::shared_ptr<BSCurveRational<T, dim>>> ||
            std::is_same_v<typename std::iterator_traits<InputIt>::value_type, std::shared_ptr<BSCurveGeneral<T, dim, false>>>, // Include this check
            "Iterator must point to elements of type BSCurve<T, dim>, BSCurveRational<T, dim>, or their shared_ptr versions"
        );
        // Determine if the curves are rational or not
        constexpr bool rational = std::is_same_v<typename std::iterator_traits<InputIt>::value_type, BSCurveRational<T, dim>>;
        // Determine we are dealing with pointers or not
        constexpr bool isPointer = std::is_pointer_v<typename std::iterator_traits<InputIt>::value_type> ||
                                    std::is_same_v<typename std::iterator_traits<InputIt>::value_type, std::shared_ptr<BSCurve<T, dim>>> ||
                                    std::is_same_v<typename std::iterator_traits<InputIt>::value_type, std::shared_ptr<BSCurveRational<T, dim>>> ||
                                    std::is_same_v<typename std::iterator_traits<InputIt>::value_type, std::shared_ptr<BSCurveGeneral<T, dim, false>>>;

        // Adjust the dimensionality based on whether the curves are rational
        std::vector<BSCurveInfo<T, dim + (rational ? 1 : 0)>> curves_info(std::distance(first, last));
        // Extract information
        auto curve_info_extractor = [](const auto& curve) {
            if constexpr (isPointer) {
                return curve->info();
            } else {
                return curve.info();
            }
        };

        std::transform(
            first, last,
            curves_info.begin(),
            curve_info_extractor
        );

        return curves_info;
    }

    /**
     * Unifies the degree of NURBS curves to the maximum degree found among them.
     * 
     * This function iterates over a collection of NURBS curves, each described by a tuple
     * containing its control points (poles), knot vector, and degree. It identifies the maximum
     * degree present among these curves and then increases the degree of each curve to match this
     * maximum degree.
     * 
     * @tparam T The floating point type of the curve coordinates and knots.
     * @tparam dim The dimensionality of the control points (poles).
     * @param curves_info A reference to a vector of tuples. Each tuple represents a NURBS curve
     *                    with its control points (poles), knot vector, and degree.
     */
    template <std::floating_point T, size_t dim>
    auto unify_degree(std::vector<gbs::BSCurveInfo<T, dim>> &curves_info) -> void
    {
        // Find the maximum degree among all curves
        auto max_degree_it = std::max_element(
            curves_info.begin(),
            curves_info.end(),
            [](const auto &C1, const auto &C2) {
                return std::get<2>(C1) < std::get<2>(C2);
            }
        );

        auto max_degree = std::get<2>(*max_degree_it);

        // Increase the degree of each curve to the maximum
        std::for_each(
            curves_info.begin(),
            curves_info.end(),
            [max_degree](auto &C) {
                auto &[poles, flat_knots, degree] = C;
                gbs::increase_degree(flat_knots, poles, degree, max_degree-degree);
                degree = max_degree;
            });
    }

    /**
     * Unifies the knots of NURBS curves by ensuring each curve has the same set of knots.
     *
     * This function iterates over the input knot-multiplicity pairs (km_in) and inserts
     * knots into the output knot vector (k_out) and control points (poles_out) as needed
     * to match the multiplicities. It also updates the output knot-multiplicity pairs (km_out).
     *
     * @tparam T Floating point type of the knots.
     * @tparam dim Dimensionality of the control points (poles).
     * @param km_out Reference to the vector of knot-multiplicity pairs for the output curve.
     * @param k_out Reference to the vector of knots for the output curve.
     * @param poles_out Reference to the vector of control points for the output curve.
     * @param degree Degree of the NURBS curve.
     * @param km_in Const reference to the vector of knot-multiplicity pairs for the input curve.
     */
    template <std::floating_point T, size_t dim>
    void unify_knots(std::vector<std::pair<T, size_t>> &km_out, std::vector<T> &k_out, std::vector<std::array<T, dim>> &poles_out, size_t degree, const std::vector<std::pair<T, size_t>> &km_in)
    {
        for (auto [u, m] : km_in) // Insert current curve knots if requierd to the first one
        {
            // Find if the knot 'u' is already present in the output knots
            auto it = std::ranges::find_if(
                km_out, 
                [u](const auto &km0_)
                { return std::abs(km0_.first - u) < gbs::knot_eps<T>; });
            if (it != km_out.end()) // Current knot is present in the first curve
            {
                if(m>it->second)
                {
                    gbs::insert_knots(u, degree, m - it->second, k_out, poles_out);
                    it->second = m;
                }
            }
            else
            {
                gbs::insert_knots(u, degree, m, k_out, poles_out);
                // Insert the knot and its multiplicity into the output knot-multiplicity pairs
                // TODO use ranges
                //    auto it_ = std::ranges::lower_bound(km_out, u, [](const auto &kmi1, auto value)
                //                                 { return kmi1.first < value; });
                auto it_ = std::lower_bound(km_out.begin(), km_out.end(), u, [](const auto &kmi1, auto value)
                                            { return kmi1.first < value; });
                km_out.insert(it_, std::make_pair(u, m));
            }
        }
    }

    /**
     * Unifies the knot vectors of a set of NURBS curves.
     *
     * This function takes a set of NURBS curves and modifies them so that they all share the same knot vector.
     * Each curve in the set is represented as a tuple containing its control points (poles), knot vector, 
     * and degree. The function first adjusts the bounds of the knot vectors of all curves to match those of
     * the first curve. It then unifies the knot vectors by ensuring each curve has the same set of knots 
     * with the same multiplicities.
     *
     * @tparam T Floating point type of the knots.
     * @tparam dim Dimensionality of the control points (poles).
     * @param curves_info A reference to a vector of tuples, where each tuple represents a NURBS curve
     *                    with its control points (std::vector<std::array<T, dim>>), 
     *                    knot vector (std::vector<T>), and degree (size_t).
     */
    template <std::floating_point T, size_t dim>
    auto unify_knots(std::vector<gbs::BSCurveInfo<T, dim>> &curves_info) -> void
    {
        // Extract the first curve's info to use as a reference
        auto &[poles0, k0, degree0] = curves_info.front();

        // Determine the start and end of the knot vector
        auto k_start = k0.front();
        auto k_end = k0.back();

        // Convert the first curve's knot vector to a multiplicity format
        auto km0 = gbs::unflat_knots(k0);

        // Adjust bounds and unify knots for each curve with the first curve
        std::for_each(
            std::next(curves_info.begin()),
            curves_info.end(),
            [k_start, k_end, degree0, &k0, &km0, &poles0](auto &C_)
            {
                auto &k = std::get<1>(C_);
                gbs::change_bounds(k_start, k_end, k);
                auto km = gbs::unflat_knots(k);
                unify_knots(km0, k0, poles0, degree0, km);
            });

        // Apply the unified knot vector to all curves
        std::for_each(
            std::next(curves_info.begin()),
            curves_info.end(),
            [&](auto &C_)
            {
                auto &[poles, k, degree] = C_;
                auto km = gbs::unflat_knots(k);
                unify_knots(km, k, poles, degree, km0);
            });
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
                insert_knots(u, degree, m, pk1.second, pk1.first);
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