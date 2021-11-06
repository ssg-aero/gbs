#pragma once
#include <gbs/basisfunctions.h>
#include <gbs/knotsfunctions.h>
#include <gbs/maths.h>
#include <gbs/vecop.h>
#include <gbs/exceptions.h>

#include <vector>
#include <array>
#include <any>
#include <type_traits>
#include <iostream>
namespace gbs
{
    template <typename T, size_t dim>
    class Geom
    {
        static_assert(std::is_floating_point<T>::value, "Only real value permitted" );
    };
    // TODO add bounded curves
    template <typename T, size_t dim>
    class Curve : public Geom<T,dim>
    {
    public:
        /**
         * @brief Curve evaluation at parameter u
         *
         * @param u : parameter on curve
         * @param d : derivative order
         * @return std::array<T, dim>
         */
        virtual auto value(T u, size_t d = 0) const -> std::array<T, dim> = 0;
        /**
         * @brief Returns curves's start stop values {U1,U2}
         * 
         * @return std::array<T,2> 
         */
        virtual auto bounds() const -> std::array<T, 2> = 0;
        auto midPosition() const -> T 
        { 
            auto [ u1, u2 ] = bounds();
            return T(0.5) * ( u1 + u2);
        }
        // virtual auto changeBounds(std::array<T, 2> &) const -> void = 0;
        // virtual auto bounded() const -> bool = 0;

        auto operator()(T u, size_t d = 0)  const -> std::array<T, dim> {return value(u,d);};
        /**
         * @brief Curve's begin
         * 
         * @param d : derivative order
         * @return std::array<T, dim> 
         */
        auto begin(size_t d = 0) const -> std::array<T, dim>
        {
            return this->value(bounds()[0], d);
        }
        /**
         * @brief Curve's end
         * 
         * @param d : derivative order
         * @return std::array<T, dim> 
         */
        auto end(size_t d = 0) const -> std::array<T, dim>
        {
            return this->value(bounds()[1], d);
        }

    };

    template<typename T>
    class OutOfBoundsCurveEval : public OutOfBoundsEval<T>
    {
        public:
        explicit OutOfBoundsCurveEval(T v, const std::array<T, 2> &bounds) : OutOfBoundsEval<T>{v, bounds, "Curve eval "} {}
    };

    /**
     * @brief check if curve's menber fullfill bspline definition
     **/ 
    template <typename T>
    auto check_curve(size_t np, const std::vector<T> &knots, size_t p)
    {
        bool ok = true;
        std::vector<size_t> m;
        std::vector<T> k;
        unflat_knots(knots,m,k);
        // needs to remove dulpicate knots, std::greater_equal<> does not support reflexivity
        // Cf. https://stackoverflow.com/questions/43696477/why-does-is-sorted-not-work-correctly-with-greater-equal
        // ok = ok && std::is_sorted(k.begin(), k.end(), std::greater<T>());
        ok = ok && std::is_sorted(k.begin(), k.end());
        if(!ok)
        {
            std::cerr << "check_curve: knots are not ordered\n";
            return ok;
        }
        ok = ok && ((p + 1 + np) == knots.size());
        if(!ok)
        {
            std::cerr << "check_curve: incorrect poles/knots/degree combination\n";
        }
        return ok;
    }
    /**
     * @brief Géneral BSpline curve class, any kind of precision, space dimension with rational definition capability
     * 
     * @tparam T        : curve precision
     * @tparam dim      : space dimension of curve (aka 1D, 2D, 3D,...)
     * @tparam rational : use last coordinate as weight
     */
    template <typename T, size_t dim, bool rational>
    class BSCurveGeneral : public Curve<T,dim>
    {
        size_t m_deg;
        std::vector<std::array<T, dim + rational>> m_poles;
        std::vector<T> m_knotsFlats;
        std::array<T,2> m_bounds;
    public:
        BSCurveGeneral() = default ;
        BSCurveGeneral( const BSCurveGeneral<T,dim,rational> &bsc ) = default ;
        /**
     * @brief Construct a new BSCurve object, non rational definition
     * 
     * @param poles : array of poles
     * @param knots : array of knots
     * @param mult  : array of knots multiplicity
     * @param deg   : curve's degree
     */
        BSCurveGeneral(const std::vector<std::array<T, dim + rational>> &poles,
                       const std::vector<T> &knots,
                       const std::vector<size_t> &mult,
                       size_t p) : m_poles(poles),
                                     m_deg(p),
                                     m_knotsFlats(flat_knots(knots, mult)),
                                     m_bounds{knots.front(),knots.back()}
        {
            if(!check_curve(poles.size(),m_knotsFlats,p)) throw std::exception("BSpline Curve constructor error.");
        }
        /**
         * @brief Construct a new BSCurve object, rational definition
         * 
         * @param poles   : array of poles
         * @param weights : array of poles' weights
         * @param knots   : array of knots
         * @param mult    : array of knots multiplicity
         * @param deg     : curve's degree
         */
        BSCurveGeneral(const std::vector<std::array<T, dim>> &poles,
                       const std::vector<T> &weights,
                       const std::vector<T> &knots,
                       const std::vector<size_t> &mult,
                       size_t p) : m_poles(merge_weights(poles, weights)),
                                     m_deg(p),
                                     m_knotsFlats(flat_knots(knots, mult)),
                                     m_bounds{knots.front(),knots.back()}
        {
            if(!check_curve(poles.size(),m_knotsFlats,p)) throw std::exception("BSpline Curve constructor error.");
        }
        /**
         * @brief Construct a new BSCurve object
         * 
         * @param poles       : array of poles
         * @param knots_flats : flat knots 
         * @param deg         : curve's degree
         */
        BSCurveGeneral(const std::vector<std::array<T, dim + rational>> &poles,
                       const std::vector<T> &knots_flats,
                       size_t p) : m_poles(poles),
                                     m_knotsFlats(knots_flats),
                                     m_deg(p),
                                     m_bounds{knots_flats.front(),knots_flats.back()}
        {
            if(!check_curve(m_poles.size(),m_knotsFlats,p)) throw std::exception("BSpline Curve constructor error.");
        }
        auto isRational() -> bool
        {
            return rational;
        }
        /**
         * @brief curve's degree
         * 
         * @return size_t 
         */
        auto degree() const noexcept -> size_t
        {
            return m_deg;
        }
        /**
         * @brief curves's flat knots
         * 
         * @return const std::vector<T>& 
         */
        auto knotsFlats() const noexcept -> const std::vector<T> &
        {
            return m_knotsFlats;
        }
        /**
         * @brief returns unflated knots
         * 
         * @return const std::vector<T>& 
         */
        auto knots() const noexcept -> const std::vector<T>
        {
            return knots_and_mults(knotsFlats()).first;
        }
        /**
         * @brief return knots multiplicities 
         * 
         * @return const std::vector<T>& 
         */
        auto mults() const noexcept -> const std::vector<size_t>
        {
            return knots_and_mults(knotsFlats()).second;
        }
        /**
         * @brief Insert knot with the given multiplicity
         * 
         * @param u : knot value
         * @param m : knot's multiplicity
         */
        auto insertKnot(T u, size_t m = 1) -> void //Fail safe, i.e. if fails, curve stays in previous state
        {
            for (auto i = 0; i < m; i++)
                insert_knot(u, m_deg, m_knotsFlats, m_poles);
        }
        /**
         * @brief Insert knots with the given multiplicities
         * 
         * @param km vector of knot/multiplicity pair
         */
        auto insertKnots(std::vector<std::pair<T,size_t>> &km) -> void
        {
            std::for_each(km.begin(), km.end(), [&](const auto km_) {
                auto m = multiplicity(m_knotsFlats, km_.first); // TODO: improve efficiency
                if (m < km_.second)
                    insertKnot(km_.first, km_.second - m);
            });
        }

        /**
         * @brief Try to remove m times the given knot
         * 
         * @param u   : knot value
         * @param tol : tolerance on curve
         * @param m   : knot removal occurrences
         */
        auto removeKnot(T u, T tol, size_t m = 1) -> void //Fail safe, i.e. if fails, curve stays in previous state
        {
            for (auto i = 0; i < m; i++)
                remove_knot(u, m_deg, m_knotsFlats, m_poles, tol);
        }
        /**
         * @brief Curve's poles
         * 
         * @return const std::vector<std::array<T, dim>>& 
         */
        auto poles() const noexcept -> const std::vector<std::array<T, dim + rational>> &
        {
            return m_poles;
        }
        constexpr auto poles_begin() const noexcept { return m_poles.begin();}
        constexpr auto poles_end()   const noexcept { return m_poles.end();}
        /**
         * @brief Copy Poles, throw std::length_error is thrown if lengths are not the same
         * 
         * @param poles 
         */
        auto copyPoles(const std::vector<std::array<T, dim + rational>> &poles) ->void
        {
            if (poles.size() != m_poles.size())
            {
                throw std::length_error("BSCurveGeneral: wrong pole vector length.");
            }
            m_poles = poles;
        }
        /**
         * @brief Move pole vector, , throw std::length_error is thrown if lengths are not the same
         * 
         * @param poles 
         */
        auto movePoles(std::vector<std::array<T, dim + rational>> &poles) ->void
        {
            if (poles.size() != m_poles.size())
            {
                throw std::length_error("BSCurveGeneral: wrong pole vector length.");
            }
            m_poles = std::move(poles);
        }
        /**
         * @brief Edit specified pole
         * 
         * @param id 
         * @param position 
         * @return auto 
         */
        auto changePole(size_t id, const point<T,dim+rational> &position)
        {
            m_poles[id] = position;
        }
        /**
         * @brief Replace curve's knots
         * 
         * @param flatKnots 
         */
        auto copyKnots(const std::vector<T> &flatKnots) -> void
        {
            if(flatKnots.size()!=m_knotsFlats.size())
                throw std::length_error("BSCurveGeneral: wrong knot vector length.");

            std::copy(flatKnots.begin(),flatKnots.end(),m_knotsFlats.begin());
            m_bounds = {flatKnots.front(),flatKnots.back()};
        }

        /**
         * @brief Reverse curve orientation
         * 
         */
        auto reverse() -> void
        {
            std::reverse(m_poles.begin(), m_poles.end());
            auto k1 = m_knotsFlats.front();
            auto k2 = m_knotsFlats.back();
            std::reverse(m_knotsFlats.begin(), m_knotsFlats.end());
            std::transform(m_knotsFlats.begin(), m_knotsFlats.end(),
                           m_knotsFlats.begin(),
                           [&](const auto k_) {
                               return k1 + k2 - k_;
                           });
        }
        /**
         * @brief Permanently trim curve between u1 and u2 by inserting knots and dropping useless ones
         * 
         * @param u1 
         * @param u2 
         */
        auto trim(T u1, T u2, bool permanently=true) -> void
        {
            m_bounds = {u1,u2};
            if (permanently)
            {
                gbs::trim(m_deg, m_knotsFlats, m_poles, u1, u2);
            }
        }
        /**
         * @brief Change parametrization to fit between k1 and k2
         * 
         * @param k1 
         * @param k2 
         */
        auto changeBounds(T k1, T k2) -> void
        {
            gbs::change_bounds(k1,k2,m_knotsFlats);
            m_bounds = {k1,k2};
        }
        /**
         * @brief Change parametrization to fit between b[0] and b[1]
         * 
         * @param b 
         */
        auto changeBounds(const std::array<T,2> &b) -> void
        {
            m_bounds = b;
            gbs::change_bounds(b[0],b[1],m_knotsFlats);
        }

        virtual auto bounds() const -> std::array<T,2> override
        {
            // return {m_knotsFlats.front(),m_knotsFlats.back()};
               return {m_bounds[0],m_bounds[1]};
        }

        auto increaseDegree() -> void
        {
            gbs::increase_degree(m_knotsFlats,m_poles,m_deg);
            m_deg++;
        }

    };

    template <typename T, size_t dim>
    class BSCurve : public BSCurveGeneral<T, dim, false>
    {
        using BSCurveGeneral<T, dim, false>::BSCurveGeneral;
    public:
        BSCurve(const BSCurveGeneral<T, dim, false> &bsc) : BSCurveGeneral<T, dim, false>(bsc.poles(), bsc.knotsFlats(), bsc.degree()) {}
        virtual auto value(T u, size_t d = 0) const -> std::array<T, dim> override
        {
            if (u < this->bounds()[0] - knot_eps|| u > this->bounds()[1] + knot_eps)
            {
                throw OutOfBoundsCurveEval(u,this->bounds());
            }
            return gbs::eval_value_simple(u, this->knotsFlats(), this->poles(), this->degree(), d);
        }

    };

    template <typename T, size_t dim>
    class BSCurveRational : public BSCurveGeneral<T, dim, true>
    {
        using  BSCurveGeneral<T, dim, true>::BSCurveGeneral;
    public:
        BSCurveRational(const BSCurveGeneral<T, dim, true> &bsc) : BSCurveGeneral<T, dim, true>(bsc.poles(), bsc.knotsFlats(), bsc.degree()) {}
        BSCurveRational(const std::vector<std::array<T, dim>> &poles,
                        const std::vector<T> &knots_flats,
                        size_t deg) : BSCurveGeneral<T, dim, true>(add_weights_coord(poles), knots_flats, deg) {}
        BSCurveRational(const BSCurve<T,dim> &crv) : BSCurveGeneral<T, dim, true>(
            add_weights_coord(crv.poles()), crv.knotsFlats(), crv.degree()) {}
        virtual auto value(T u, size_t d = 0) const -> std::array<T, dim> override
        {            
            if (u < this->bounds()[0] - knot_eps|| u > this->bounds()[1] + knot_eps)
            {
                throw OutOfBoundsCurveEval(u,this->bounds());
            }
            return eval_rational_value_simple<T,dim>(u,this->knotsFlats(),this->poles(),this->degree(),d);
        }
        auto polesProjected() const -> points_vector<T,dim>
        {
            points_vector<T,dim> poles_(this->poles().size());
            std::transform(this->poles().begin(),this->poles().end(),poles_.begin(),
            [](const auto &p){return weight_projection(p);});
            return poles_;
        }

        auto weights() const -> std::vector<T>
        {
            std::vector<T> weights_(this->poles().size());
            std::transform(this->poles().begin(),this->poles().end(),weights_.begin(),
            [](const auto &p){return p.back();});
            return weights_;
        }

    };

    template< typename T>
    class BSCfunction
    {
        BSCurve<T,1> crv_;
        public:
        BSCfunction(const BSCfunction<T> &func) : crv_{func.basisCurve()} 
        {}
        BSCfunction(const BSCurve<T,1> &crv) : crv_{crv} 
        {}
        BSCfunction(const std::vector<std::array<T, 1>> &poles, const std::vector<T> &knots_flats, size_t deg) : crv_{poles,knots_flats,deg}
        {}
        BSCfunction(const std::vector<T> &poles, const std::vector<T> &knots_flats, size_t deg) : crv_{poles_array(poles),knots_flats,deg}
        {}
        BSCfunction(const std::vector<T> &poles, const std::vector<T> &knots, const std::vector<size_t> mults, size_t deg) : crv_{poles_array(poles), knots, mults, deg}
        {}
        BSCfunction() = default;
        auto operator()(T u, size_t d = 0)  const -> T
        {
            return value(u,d);
        }
        auto value(T u, size_t d = 0)  const -> T
        {
            return crv_.value(u,d)[0];
        }
        auto bounds() const ->  std::array<T, 2> {return crv_.bounds();}
        /**
         * @brief Change parametrization to fit between k1 and k2
         * 
         * @param k1 
         * @param k2 
         */
        auto changeBounds(T k1, T k2) -> void { crv_.changeBounds(k1,k2);}
        /**
         * @brief Change parametrization to fit between b[0] and b[1]
         * 
         * @param b 
         */
        auto changeBounds(const std::array<T,2> &b) -> void {crv_.changeBounds(b);}
        auto begin() const ->  T {return crv_.begin()[0];}
        auto end() const ->  T {return crv_.end()[0];}
        auto basisCurve() const -> const BSCurve<T,1> & {return crv_;}
        static auto poles_array(const std::vector<T> &poles) -> std::vector<std::array<T, 1>>
        {
            std::vector<std::array<T, 1>> arr(poles.size());
            std::transform(
                poles.begin(),
                poles.end(),
                arr.begin(),
                [](const auto &pole){return std::array<T, 1>{pole};}
            );
            return arr;
        }
    };

    template <typename T, size_t dim>



    using BSCurve2d_f = BSCurve<float,2>;
    using BSCurve3d_f = BSCurve<float,3>;
    using BSCurve2d_d = BSCurve<double,2>;
    using BSCurve3d_d = BSCurve<double,3>;
    using BSCurveRational2d_f = BSCurveRational<float,2>;
    using BSCurveRational3d_f = BSCurveRational<float,3>;
    using BSCurveRational2d_d = BSCurveRational<double,2>;
    using BSCurveRational3d_d = BSCurveRational<double,3>;
} // namespace gbs