#pragma once
#include <Geom2d_BSplineCurve.hxx>
#include <Geom_BSplineCurve.hxx>
#include <TColgp_Array1OfPnt2d.hxx>
#include <TColgp_Array1OfPnt.hxx>
#include <TColgp_Array1OfVec2d.hxx>
#include <TColgp_Array1OfVec.hxx>
#include <vector>
#include <occt-utils/containers.h>

#include <gbslib/bscurve.h>
#include <gbslib/knotsfunctions.h>
#include <tuple>
namespace occt_utils
{
    template <typename T, size_t dim>
    auto extract_weights(const std::vector<std::array<T, dim>> &poles) -> std::tuple<std::vector<std::array<T, dim - 1>>, std::vector<T>>
    {
        std::vector<std::array<T, dim - 1>> p(poles.size());
        std::vector<T> w(p.size());

        std::transform(
            std::execution::par,
            poles.begin(), poles.end(),
            w.begin(),
            [&](const auto &p_) { return p_.back(); });

        std::transform(
            std::execution::par,
            poles.begin(), poles.end(),
            w.begin(),
            p.begin(),
            [&](const auto &p_, const auto &w_) {
                std::array<T, dim - 1> pole;
                std::transform(
                    p_.begin(), std::next(p_.end(), -1), pole.begin(),
                    [&](const auto coord_) { return coord_ / w_; });
                return pole;
            });

        return std::make_tuple(p, w);
    }

    auto inline BSplineCurve(  const std::vector<gp_Pnt2d> &p,
                        const std::vector<double> &k,
                        const std::vector<int> &m,
                        int deg) -> Handle(Geom2d_BSplineCurve)
    { return new Geom2d_BSplineCurve(
        occt_utils::col_from_vec(p),
        occt_utils::col_from_vec(k),
        occt_utils::col_from_vec(m),deg);}

    auto inline BSplineCurve(  const std::vector<gp_Pnt> &p,
                        const std::vector<double> &k,
                        const std::vector<int> &m,
                        int deg) -> Handle(Geom_BSplineCurve)
    { return new Geom_BSplineCurve(
        occt_utils::col_from_vec(p),
        occt_utils::col_from_vec(k),
        occt_utils::col_from_vec(m),deg);}

    auto inline BSplineCurve(  const std::vector<std::array<double,2> > &p,
                        const std::vector<double> &k,
                        const std::vector<int> &m,
                        int deg) -> Handle(Geom2d_BSplineCurve)
    { return new Geom2d_BSplineCurve(
        occt_utils::col_from_vec<gp_Pnt2d>(p),
        occt_utils::col_from_vec(k),
        occt_utils::col_from_vec(m),deg);}

    auto inline BSplineCurve(  const std::vector<std::array<double,3> > &p,
                        const std::vector<double> &k,
                        const std::vector<int> &m,
                        int deg) -> Handle(Geom_BSplineCurve)
    { return new Geom_BSplineCurve(
        occt_utils::col_from_vec<gp_Pnt>(p),
        occt_utils::col_from_vec(k),
        occt_utils::col_from_vec(m),deg);}

    auto inline BSplineCurve(  const std::vector<std::array<double,2> > &p,
                        const std::vector<double> &w,
                        const std::vector<double> &k,
                        const std::vector<int> &m,
                        int deg) -> Handle(Geom2d_BSplineCurve)
    { return new Geom2d_BSplineCurve(
        occt_utils::col_from_vec<gp_Pnt2d>(p),
        occt_utils::col_from_vec<double>(w),
        occt_utils::col_from_vec(k),
        occt_utils::col_from_vec(m),deg);}

    auto inline BSplineCurve(  const std::vector<std::array<double,3> > &p,
                        const std::vector<double> &w,
                        const std::vector<double> &k,
                        const std::vector<int> &m,
                        int deg) -> Handle(Geom_BSplineCurve)
    { return new Geom_BSplineCurve(
        occt_utils::col_from_vec<gp_Pnt>(p),
        occt_utils::col_from_vec<double>(w),
        occt_utils::col_from_vec(k),
        occt_utils::col_from_vec(m),deg);}

    auto inline BSplineCurve(  const std::vector<std::array<double,4> > &p,
                        const std::vector<double> &k,
                        const std::vector<int> &m,
                        int deg) -> Handle(Geom_BSplineCurve)
    { 

        auto poles_and_weigts = extract_weights(p);

        return new Geom_BSplineCurve(
        occt_utils::col_from_vec<gp_Pnt>(std::get<0>(poles_and_weigts)),
        occt_utils::col_from_vec<double>(std::get<1>(poles_and_weigts)),
        occt_utils::col_from_vec(k),
        occt_utils::col_from_vec(m),deg);}

    auto inline BSplineCurve( const gbs::BSCurve<double,2> &c ) -> Handle(Geom2d_BSplineCurve)
    {
        std::vector<int> mult;
        std::vector<double> knots;
        gbs::unflat_knots(c.knotsFlats(), mult, knots);
        return occt_utils::BSplineCurve(c.poles(), knots, mult, c.degree());
    }

    auto inline NURBSplineCurve( const gbs::BSCurve<double,3> &c ) -> Handle(Geom2d_BSplineCurve)
    {
        std::vector<int> mult;
        std::vector<double> knots;
        gbs::unflat_knots(c.knotsFlats(), mult, knots);
        auto poles_and_weigts = extract_weights(c.poles());
        auto poles   = std::get<0>(poles_and_weigts);
        auto weights = std::get<1>(poles_and_weigts);
        return occt_utils::BSplineCurve(poles, weights, knots, mult, c.degree());
    }

    auto inline BSplineCurve( const gbs::BSCurve<double,3> &c ) -> Handle(Geom_BSplineCurve)
    {
        std::vector<int> mult;
        std::vector<double> knots;
        gbs::unflat_knots(c.knotsFlats(), mult, knots);
        return occt_utils::BSplineCurve(c.poles(), knots, mult, c.degree());
    }

    auto inline NURBSplineCurve( const gbs::BSCurve<double,4> &c ) -> Handle(Geom_BSplineCurve)
    {
        std::vector<int> mult;
        std::vector<double> knots;
        gbs::unflat_knots(c.knotsFlats(), mult, knots);
        auto poles_and_weigts = extract_weights(c.poles());
        auto poles   = std::get<0>(poles_and_weigts);
        auto weights = std::get<1>(poles_and_weigts);
        return occt_utils::BSplineCurve(poles, weights, knots, mult, c.degree());
    }

    Standard_EXPORT auto to_bs_list(const std::list<Handle(Geom_Curve)> &lst) -> std::list<Handle(Geom_Curve)>;

    Standard_EXPORT auto bs_from_list(const std::list<Handle(Geom_Curve)> &lst,double tol) ->  Handle(Geom_BSplineCurve);

    Standard_EXPORT auto bscurve_c1(    const TColgp_Array1OfPnt2d &pt_lst, 
                                        const gp_Vec2d &t1,
                                        const gp_Vec2d &t2,
                                        bool scale_tg,
                                        double tol
                                        ) -> Handle(Geom2d_BSplineCurve);
    Standard_EXPORT auto bscurve_c1( const TColgp_Array1OfPnt &pt_lst,
                                     const gp_Vec &t1,
                                     const gp_Vec &t2,
                                     bool scale_tg,
                                     double tol
                                     ) -> Handle(Geom_BSplineCurve);
    template<size_t d>
    auto bscurve_c1(    const std::vector< std::array<double,d> > &pt_lst, 
                        const std::array<double,d> &t1,
                        const std::array<double,d> &t2,
                        bool scale_tg,
                        double tol)
    {
        return bscurve_c1(
            col_from_vec< std::conditional<d == 2,gp_Pnt2d,gp_Pnt>::type >(pt_lst),
            vector(t1),
            vector(t2),
            scale_tg,
            tol
        );
    }

    Standard_EXPORT auto bscurve_c1(    const TColgp_Array1OfPnt2d &pt_lst, 
                                        const TColgp_Array1OfVec2d &tg_lst,
                                        bool scale_tg,double tol
                                        ) -> Handle(Geom2d_BSplineCurve);
    Standard_EXPORT auto bscurve_c1(    const TColgp_Array1OfPnt &pt_lst,
                                        const TColgp_Array1OfVec &tg_lst,
                                        bool scale_tg,
                                        double tol
                                        ) -> Handle(Geom_BSplineCurve);
    template<size_t d>
    auto bscurve_c1(    const std::vector< std::array<double,d> > &pt_lst, 
                        const std::vector< std::array<double,d> > &tg_lst,
                        bool scale_tg,
                        double tol)
    {
        return bscurve_c1(
            col_from_vec< std::conditional<d == 2,gp_Pnt2d,gp_Pnt>::type >(pt_lst),
            col_from_vec< std::conditional<d == 2,gp_Vec2d,gp_Vec>::type >(tg_lst),
            scale_tg,
            tol
        );
    }
    Standard_EXPORT auto bscurve_c1(    const std::vector< std::array<double,2> > &pt_lst, 
                                        const std::vector< double > &pr_lst, 
                                        const std::array<double,2> &t1,
                                        const std::array<double,2> &t2,
                                        bool scale_tg,double tol
                                        ) -> Handle(Geom2d_BSplineCurve);
    Standard_EXPORT auto bscurve_c1(    const std::vector< std::array<double,2> > &pt_lst,
                                        const std::vector< double > &pr_lst,
                                        const std::array<double,2> &t1,
                                        const std::array<double,2> &t2,
                                        bool scale_tg,
                                        double tol
                                        ) -> Handle(Geom2d_BSplineCurve);
    Standard_EXPORT auto bscurve_c1(    const std::vector< std::array<double,2> > &pt_lst,
                                        double tol
                                        ) -> Handle(Geom2d_BSplineCurve);
    Standard_EXPORT auto bscurve_c1(    const Handle(Geom2d_Curve) &back,
                                        const Handle(Geom2d_Curve) &front,
                                        bool scale_tg,
                                        double tol
                                        ) -> Handle(Geom2d_BSplineCurve);

    Standard_EXPORT auto bscurve_c1(    const TColgp_Array1OfPnt &pt_lst,
                                        const TColgp_Array1OfVec &tg_lst,
                                        bool scale_tg,
                                        double tol
                                        ) -> Handle(Geom_BSplineCurve);
    
    Standard_EXPORT auto bscurve_c2_approx( const TColgp_Array1OfPnt2d &pt_lst,
                                            const TColgp_Array1OfVec2d &tg_lst,
                                            double tol
                                            ) -> Handle(Geom2d_BSplineCurve);
    Standard_EXPORT auto bscurve_c2_approx( const TColgp_Array1OfPnt &pt_lst,
                                            const TColgp_Array1OfVec &tg_lst,
                                            double tol
                                            ) -> Handle(Geom_BSplineCurve);
    template<size_t d>
    auto bscurve_c2_approx( const std::vector< std::array<double,d> > &pt_lst,
                            const std::vector< std::array<double,d> > &tg_lst,
                            double tol
                            )
    {
        return bscurve_c2_approx(
            col_from_vec< std::conditional<d == 2,gp_Pnt2d,gp_Pnt>::type >(pt_lst),
            col_from_vec< std::conditional<d == 2,gp_Vec2d,gp_Vec>::type >(tg_lst),
            tol
        );
    }
} // namespace occt_utils