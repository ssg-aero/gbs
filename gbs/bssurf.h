#pragma once
import basis_functions;
import knots_functions;
#include <gbs/bscurve.h>
#include <gbs/polestools.h>
#include <gbs/exceptions.h>
#include <vector>
#include <array>
#include <algorithm>
#include <ranges>
namespace gbs
{
    // TODO: move to pole surf utils file, as well as transpose_poles form loftBase.h
    /**
     * @brief 
     * 
     * @tparam T 
     * @tparam dim 
     * @param poles 
     * @param n_new_poles_u : number of new u poles
     * @return points_vector<T, dim> 
     */
    template <typename T, size_t dim>
    auto inverted_uv_poles(const points_vector<T, dim> &poles, size_t n_new_poles_u) -> points_vector<T, dim>
    {
        auto n_poles_v = poles.size() / n_new_poles_u;
        points_vector<T, dim> poles_t(poles.size());
        for (int i = 0; i < n_new_poles_u; i++)
        {
            for (int j = 0; j < n_poles_v; j++)
            {
                poles_t[i + n_new_poles_u * j] = poles[j + n_poles_v * i];
            }
        }
        return poles_t;
    }
    /**
     * @brief 
     * 
     * @tparam T 
     * @tparam dim 
     * @param poles 
     * @param n_new_poles_u : number of new u poles
     */
    template <typename T, size_t dim>
    auto invert_uv_poles(points_vector<T, dim> &poles, size_t n_new_poles_u) -> void
    {
        poles = std::move(inverted_uv_poles(poles, n_new_poles_u));
    }

    template<typename T>
    class OutOfBoundsSurfaceUEval : public OutOfBoundsEval<T>
    {
        public:
        explicit OutOfBoundsSurfaceUEval(T v, const std::array<T, 2> &bounds) : OutOfBoundsEval<T>{v, bounds, "Surface U eval "} {}
    };

    template<typename T>
    class OutOfBoundsSurfaceVEval : public OutOfBoundsEval<T>
    {
        public:
        explicit OutOfBoundsSurfaceVEval(T v, const std::array<T, 2> &bounds) : OutOfBoundsEval<T>{v, bounds, "Surface V eval "} {}
    };

    template<typename T>
    class OutOfBoundsSurfaceIsoU : public OutOfBoundsEval<T>
    {
        public:
        explicit OutOfBoundsSurfaceIsoU(T u, const std::array<T, 2> &bounds) : OutOfBoundsEval<T>{u, bounds, "Surface U iso "} {}
    };

    template<typename T>
    class OutOfBoundsSurfaceIsoV : public OutOfBoundsEval<T>
    {
        public:
        explicit OutOfBoundsSurfaceIsoV(T v, const std::array<T, 2> &bounds) : OutOfBoundsEval<T>{v, bounds, "Surface V iso "} {}
    };

    /**
     * @brief The Surface class is an abstract class representing a parametric surface.
     *
     * @tparam T The type of the coordinates.
     * @tparam dim The dimension of the space in which the surface is defined.
     */
    template <typename T, size_t dim>
    class Surface : public Geom<T, dim>
    {
    public:
        /**
         * @brief  Surface evaluation at parameters {u,v}
         *
         * @param u  The u parameter on surface.
         * @param v  The v parameter on surface.
         * @param du The u derivative order.
         * @param dv The v derivative order.
         * @return A const reference to the resulting point<T, dim>.
         */
        virtual auto value(T u, T v, size_t du = 0, size_t dv = 0) const -> point<T, dim> = 0;

        /**
         * @brief  Surface evaluation at parameters {u,v}.
         *
         * @param uv The {u,v} parameter on surface.
         * @param du The u derivative order.
         * @param dv The v derivative order.
         * @return A const reference to the resulting point<T, dim>.
         */
        auto value(const std::array<T, 2> &uv, size_t du = 0, size_t dv = 0) const -> point<T, dim>
        {
            auto [u, v] = uv;
            return value(u, v, du, dv);
        };
        /**
         * @brief Evaluates the surface at multiple parameter pairs (u, v).
         *
         * This method evaluates the surface at each pair of parameters specified in the container `uv_lst`.
         * Each element of `uv_lst` is a `std::array<T, 2>`, representing a pair of parameters (u, v).
         * The `value` function is applied to each element of `uv_lst` to compute the surface's
         * value (and optionally derivatives) at these parameter pairs. The results are stored in a
         * `std::vector<std::array<T, dim>>`.
         *
         * The template parameter `_Container` represents any container type that holds elements of type `std::array<T, 2>`.
         * The container must support `begin()` and `end()` methods.
         *
         * @tparam _Container A template template parameter representing the container type.
         * @tparam Args Additional template arguments for the container (e.g., custom allocators).
         * @param uv_lst A container of parameter pairs (u, v) at which to evaluate the surface.
         * @param du The derivative order with respect to u. Default is 0.
         * @param dv The derivative order with respect to v. Default is 0.
         * @return std::vector<std::array<T, dim>> A container with the surface evaluations (and derivatives if du > 0 or dv > 0).
         */
        template <template<typename ...> class _Container, typename... Args>
        auto values(const _Container<std::array<T, 2>, Args...> &uv_lst, size_t du = 0, size_t dv = 0) const -> points_vector<T, dim>
        {
            std::vector<std::array<T, dim>> values(uv_lst.size());
            // Parallelization doesn't seems to be worth
            std::transform(
                uv_lst.begin(), uv_lst.end(),
                values.begin(),
                [this, du, dv](const auto &uv){return this->value(uv, du, dv);}
            );
            return values;
        }
        /**
         * @brief Evaluates the surface at multiple parameter pairs (u, v) provided in separate containers.
         *
         * This method evaluates the surface at parameter pairs formed by combining elements from `u_lst` and `v_lst`.
         * The `value` function is applied to each pair of parameters (u from `u_lst`, v from `v_lst`) to compute
         * the surface's value (and optionally derivatives) at these parameter pairs. The results are stored in a
         * `std::vector<std::array<T, dim>>`.
         *
         * The template parameter `_Container` represents any container type that holds elements of type `T`.
         * Both containers must be of the same size and support `begin()` and `end()` methods.
         *
         * @tparam _Container A template template parameter representing the container type.
         * @tparam Args Additional template arguments for the container (e.g., custom allocators).
         * @param u_lst A container of u parameters.
         * @param v_lst A container of v parameters. Must be the same size as `u_lst`.
         * @param du The derivative order with respect to u. Default is 0.
         * @param dv The derivative order with respect to v. Default is 0.
         * @return std::vector<std::array<T, dim>> A container with the surface evaluations (and derivatives if du > 0 or dv > 0).
         */
        template <template<typename ...> class _Container, typename... Args>
        auto values(const _Container<T, Args...>& u_lst, const _Container<T, Args...>& v_lst, size_t du = 0, size_t dv = 0) const -> std::vector<std::array<T, dim>>
        {
            assert( u_lst.size() == v_lst.size() );
            std::vector<std::array<T, dim>> values(u_lst.size());
            // Parallelization doesn't seems to be worth
            std::transform(
                u_lst.begin(), u_lst.end(),
                v_lst.begin(),
                values.begin(),
                [this, du, dv](T u, T v){return this->value(u, v, du, dv);}
            );
            return values;
        }

        /**
         * @brief  Returns the surface's bounds {U1,U2,V1,V2}.
         *
         * @return An array containing the surface's bounds.
         */
        virtual auto bounds() const -> std::array<T, 4> = 0;

        /**
         * @brief Returns the surface's U bounds.
         *
         * @return An array containing the surface's U bounds.
         */
        auto boundsU() const -> std::array<T, 2> { return {bounds()[0], bounds()[1]}; }

        /**
         * @brief Returns the surface's V bounds.
         *
         * @return An array containing the surface's V bounds.
         */
        auto boundsV() const -> std::array<T, 2> { return {bounds()[2], bounds()[3]}; }

        /**
         * @brief Returns the surface's value at parameters {u,v}.
         *
         * @param u  The u parameter on surface.
         * @param v  The v parameter on surface.
         * @param du The u derivative order.
         * @param dv The v derivative order.
         * @return A const reference to the resulting point<T, dim>.
         */
        auto operator()(T u, T v, size_t du = 0, size_t dv = 0) const -> point<T, dim> { return value(u, v, du, dv); };

        /**
         * @brief Returns the surface's value at parameters {u,v}.
         *
         * @param uv The {u,v} parameter on surface.
         * @param du The u derivative order.
         * @param dv The v derivative order.
         * @return A const reference to the resulting point<T, dim>.
         */
        auto operator()(const std::array<T, 2> &uv, size_t du = 0, size_t dv = 0) const -> point<T, dim> { return value(uv, du, dv); };

        /**
         * @brief Computes the tangent vector in the direction of the given (delta_u, delta_v) values at the point (u, v)
         * 
         * @param u         The u parameter of the point
         * @param v         The v parameter of the point
         * @param delta_u   The delta_u parameter for the tangent direction
         * @param delta_v   The delta_v parameter for the tangent direction
         * 
         * @return The tangent vector at the given point and direction
         */
        auto tangent(T u, T v, T delta_u, T delta_v) const -> point<T, dim>
        {
            return value(u, v, 1, 0) * delta_u + value(u, v, 0, 1) * delta_v;
        }
    };

    template <std::floating_point T, size_t dim>
    using BSSurfaceInfo = std::tuple< std::vector< std::vector< std::array<T, dim> > >, std::vector<T>, std::vector<T>, size_t, size_t>;

// Forward declaration needed for isoU/V methods
    template <typename T, size_t dim>
    class BSSurface;

    template <typename T, size_t dim>
    class BSSurfaceRational;

    /**
     * @brief Géneral BSpline surface class, any kind of precision, space dimension with rational definition capability
     * 
     * @tparam T 
     * @tparam dim 
     */
    template <typename T, size_t dim, bool rational>
    class BSSurfaceGeneral : public Surface<T, dim>
    {
        size_t m_degU, m_degV;
        std::vector<std::array<T, dim + rational>> m_poles;
        std::vector<T> m_knotsFlatsU;
        std::vector<T> m_knotsFlatsV;
        std::array<T, 4> m_bounds;

    public:
        /**
     * @brief Construct a new BSSurface object, non rational definition
     * 
     * @param poles 
     * @param knots_flatsU 
     * @param knots_flatsV 
     * @param degU 
     * @param degV 
     */
        BSSurfaceGeneral(const std::vector<std::array<T, dim + rational>> &poles,
                         const std::vector<T> &knots_flatsU,
                         const std::vector<T> &knots_flatsV,
                         size_t degU,
                         size_t degV) : m_poles(poles),
                                        m_knotsFlatsU(knots_flatsU),
                                        m_knotsFlatsV(knots_flatsV),
                                        m_degU(degU),
                                        m_degV(degV),
                                        m_bounds{knots_flatsU.front(), knots_flatsU.back(), knots_flatsV.front(), knots_flatsV.back()}
        //,
        //   m_rational(false)

        {
            if (nPolesU()*nPolesV()!=m_poles.size())
                throw std::invalid_argument("BSpline Surface constructor error. Pole size error.");
            if (!check_curve(nPolesU(), m_knotsFlatsU, degU))
                throw std::invalid_argument("BSpline Surface constructor error. Invalid knotsU/degreeU combination");
            if (!check_curve(nPolesV(), m_knotsFlatsV, degV))
                throw std::invalid_argument("BSpline Surface constructor error. Invalid knotsV/degreeV combination");
        }

        BSSurfaceGeneral(const std::vector<std::array<T, dim + rational>> &poles,
                         const std::vector<T> &knotsU,
                         const std::vector<T> &knotsV,
                         const std::vector<size_t> &multsU,
                         const std::vector<size_t> &multsV,
                         size_t degU,
                         size_t degV) : BSSurfaceGeneral{poles, flat_knots(knotsU, multsU), flat_knots(knotsV, multsV), degU, degV} {}
        BSSurfaceGeneral(const BSSurfaceInfo<T, dim> &srf_info) : BSSurfaceGeneral{
                                                                      flatten_poles(std::get<0>(srf_info)),
                                                                      std::get<1>(srf_info),
                                                                      std::get<2>(srf_info),
                                                                      std::get<3>(srf_info),
                                                                      std::get<4>(srf_info)}
        { }
        BSSurfaceGeneral() = default;
        BSSurfaceGeneral(const BSSurfaceGeneral<T, dim, rational> &) = default;
        // BSSurfaceGeneral<T, dim, rational> &operator=(BSSurfaceGeneral<T, dim, rational> &srf) const = default;

        /**
         * @brief Retrieves comprehensive information about the BSpline surface.
         * 
         * This method constructs and returns a tuple containing detailed information about the BSpline surface,
         * including the control points (poles), knot vectors for both U and V directions, and the degrees in both
         * the U and V directions.
         * 
         * @return BSSurfaceInfo<T, dim+rational> A tuple containing the control points, knot vectors, and degrees.
         */
        auto info() const -> BSSurfaceInfo<T, dim+rational>
        {
            std::vector< std::vector< std::array<T, dim> > > polesUV(this->nPolesV());
            for(size_t j{}; j < this->nPolesV(); j++)
            {
                polesUV[j] = this->polesU(j);
            }
            return std::make_tuple(
                polesUV, this->knotsFlatsU(), this->knotsFlatsV(), this->degreeU(), this->degreeV()
            );
        }

        /**
         * @brief Gets the degree of the surface in the U direction.
         * 
         * @return size_t The degree of the BSpline surface in the U direction.
         */
        auto degreeU() const -> size_t
        {
            return m_degU;
        }

        /**
         * @brief Gets the degree of the surface in the V direction.
         * 
         * @return size_t The degree of the BSpline surface in the V direction.
         */
        auto degreeV() const -> size_t
        {
            return m_degV;
        }

        /**
         * @brief Gets the order of the surface in the U direction.
         * 
         * The order is defined as the degree plus one.
         * 
         * @return size_t The order of the BSpline surface in the U direction.
         */
        auto orderU() const -> size_t
        {
            return m_degU+1;
        }

        /**
         * @brief Gets the order of the surface in the V direction.
         * 
         * The order is defined as the degree plus one.
         * 
         * @return size_t The order of the BSpline surface in the V direction.
         */
        auto orderV() const -> size_t
        {
            return m_degV+1;
        }

        /**
         * @brief Gets the flattened knot vector for the U direction.
         * 
         * @return const std::vector<T> & A reference to the flattened knot vector in the U direction.
         */
        auto knotsFlatsU() const -> const std::vector<T> &
        {
            return m_knotsFlatsU;
        }

        /**
         * @brief Gets the flattened knot vector for the V direction.
         * 
         * This method returns a reference to the internal knot vector in the V direction. The knot vector
         * is "flattened", meaning it's presented as a single-dimensional array.
         * 
         * @return const std::vector<T> & A reference to the flattened knot vector in the V direction.
         */
        auto knotsFlatsV() const -> const std::vector<T> &
        {
            return m_knotsFlatsV;
        }
        /**
         * @brief returns unflated knots U
         * 
         * @return const std::vector<T>& 
         */
        auto knotsU() const -> const std::vector<T>
        {
            return knots_and_mults(knotsFlatsU()).first;
        }
        /**
         * @brief returns unflated knots V
         * 
         * @return const std::vector<T>& 
         */
        auto knotsV() const -> const std::vector<T>
        {
            return knots_and_mults(knotsFlatsV()).first;
        }
        /**
         * @brief return knots U multiplicities 
         * 
         * @return const std::vector<T>& 
         */
        auto multsU() const -> const std::vector<size_t>
        {
            return knots_and_mults(knotsFlatsU()).second;
        }
        /**
         * @brief return knots V multiplicities 
         * 
         * @return const std::vector<T>& 
         */
        auto multsV() const -> const std::vector<size_t>
        {
            return knots_and_mults(knotsFlatsV()).second;
        }

        /**
         * @brief Inserts a knot into the U direction of the BSpline surface.
         * 
         * This method inserts a knot 'u' into the surface, repeating 'mult' times. 
         * It updates the control points and knot vector accordingly. The operation is fail-safe, 
         * meaning if it fails, the surface remains in its previous state.
         * 
         * @param u The knot value to insert.
         * @param mult The multiplicity of the knot to insert.
         * @return size_t The index at which the knot was inserted.
         */
        auto insertKnotU(T u, size_t mult = 1) //Fail safe, i.e. if fails, surf stays in previous state
        {
            auto m = nPolesV();

            points_vector<T, dim + rational> poles_;
            std::vector<T> knots_flatsU_;

            size_t ik{};
            for (auto j = 0; j < m; j++)
            {
                auto polesU_ = polesU(j);
                knots_flatsU_ = m_knotsFlatsU;

                ik = insert_knots(u, m_degU, mult, knots_flatsU_, polesU_);

                poles_.insert(poles_.end(), polesU_.begin(), polesU_.end());
            }

            m_knotsFlatsU = std::move(knots_flatsU_);
            m_poles = std::move(poles_);
            return ik;
        }

        /**
         * @brief Inserts a knot into the V direction of the BSpline surface.
         * 
         * Similar to insertKnotU, this method inserts a knot 'u' into the V direction of the surface, 
         * repeating 'mult' times, and updates the control points and knot vector. The operation 
         * is also fail-safe.
         * 
         * @param u The knot value to insert.
         * @param mult The multiplicity of the knot to insert.
         * @return size_t The index at which the knot was inserted.
         */
        auto insertKnotV(T u, size_t mult = 1) //Fail safe, i.e. if fails, surf stays in previous state
        {
            auto n = nPolesU();

            points_vector<T, dim + rational> poles_;
            std::vector<T> knots_flatsV_;

            size_t ik{};
            for (auto j = 0; j < n; j++)
            {
                auto polesV_ = polesV(j);
                knots_flatsV_ = m_knotsFlatsV;
                ik = insert_knots(u, m_degV, mult, knots_flatsV_, polesV_);

                poles_.insert(poles_.end(), polesV_.begin(), polesV_.end());
            }

            invert_uv_poles(poles_, nPolesU());

            m_knotsFlatsV = std::move(knots_flatsV_);
            m_poles = std::move(poles_);
            return ik;
        }

        /**
         * @brief Provides access to the control points of the surface.
         * 
         * @return const points_vector<T, dim + rational>& A reference to the control points.
         */
        auto poles() const noexcept -> const points_vector<T, dim + rational> &
        {
            return m_poles;
        }

        /**
         * @brief Accesses a specific control point on the surface by its indices.
         * 
         * @param i The index in the U direction.
         * @param j The index in the V direction.
         * @return point<T, dim + rational>& A reference to the specified control point.
         * @throws std::out_of_range If the indices are out of bounds.
         */
        auto pole(size_t i, size_t j) -> point<T, dim + rational> &
        {
            auto n = nPolesU();
            auto id = j + i * n;
            if(id>=m_poles.size())
            {
                throw std::out_of_range("BSSurfaceGeneral::pole(size_t i, size_t j) out of bounds index.");
            }
            return m_poles[id];
        }

        /**
         * @brief Const overload of the pole access method.
         * 
         * @param i The index in the U direction.
         * @param j The index in the V direction.
         * @return const point<T, dim + rational>& A const reference to the specified control point.
         * @throws std::out_of_range If the indices are out of bounds.
         */
        auto pole(size_t i, size_t j) const -> const point<T, dim + rational> &
        {
            auto n = nPolesU();
            auto id = j + i * n;
            if(id>=m_poles.size())
            {
                throw std::out_of_range("BSSurfaceGeneral::pole(size_t i, size_t j) out of bounds index.");
            }
            return m_poles[id];
        }

        /**
         * @brief Retrieves the control points in the U direction for a given V index.
         * 
         * This method constructs and returns a vector of control points in the U direction for a specified index in the V direction.
         * The method is guaranteed not to throw exceptions. Last PolesU vector is sent if overflow.
         * 
         * @param j The index in the V direction.
         * @return const points_vector<T, dim + rational> A vector of control points in the U direction.
         */
        auto polesU(size_t j) const noexcept -> const points_vector<T, dim + rational>
        {
            auto n = nPolesU();
            auto m = nPolesV();
            points_vector<T, dim + rational> poles_(n);
            auto beg_ = std::next(m_poles.begin(), j * n);
            auto end_ = (j < m) ? std::next(m_poles.begin(), (j + 1) * n) : m_poles.end();
            std::copy(beg_, end_, poles_.begin());

            return poles_;
        }

        /**
         * @brief Retrieves the control points in the V direction for a given U index.
         * 
         * This method constructs and returns a vector of control points in the V direction for a specified index in the U direction.
         * It iterates over the control points array, selecting appropriate points for the specified V direction.
         * The method is guaranteed not to throw exceptions. Last PolesV vector is sent if overflow.
         * 
         * @param j The index in the U direction.
         * @return const points_vector<T, dim + rational> A vector of control points in the V direction.
         */
        auto polesV(size_t i) const noexcept -> const points_vector<T, dim + rational>
        {
            auto n = nPolesU();
            auto m = nPolesV();
            points_vector<T, dim + rational> poles_(m);
            for (auto j = 0; j < m; j++)
            {
                poles_[j] = m_poles[i + j * n];
            }

            return poles_;
        }

        /**
         * @brief Calculates the number of control points in the U direction.
         * 
         * @return size_t The number of control points in the U direction.
         */
        auto nPolesU() const -> size_t
        {
            return m_knotsFlatsU.size() - m_degU - 1;
        }

        /**
         * @brief Calculates the number of control points in the V direction.
         * 
         * @return size_t The number of control points in the V direction.
         */
        auto nPolesV() const -> size_t
        {
            return m_knotsFlatsV.size() - m_degV - 1;
        }

        /**
         * @brief Provides an iterator to the beginning of the control points.
         * 
         * @return An iterator to the start of the control points.
         */
        constexpr auto poles_begin() const noexcept { return m_poles.begin(); }

        /**
         * @brief Provides an iterator to the end of the control points.
         * 
         * @return An iterator to the end of the control points.
         */
        constexpr auto poles_end() const noexcept { return m_poles.end(); }
        /**
         * @brief Copy Poles, throw std::length_error is thrown if lengths are not the same
         * 
         * @param poles 
         */
        auto copyPoles(const std::vector<std::array<T, dim + rational>> &poles) -> void
        {
            if (poles.size() != m_poles.size())
            {
                throw std::length_error("BSSurfaceGeneral: wrong pole vector length.");
            }
            m_poles = poles;
        }
        /**
         * @brief Move pole vector, , throw std::length_error is thrown if lengths are not the same
         * 
         * @param poles 
         */
        auto movePoles(std::vector<std::array<T, dim + rational>> &poles) -> void
        {
            if (poles.size() != m_poles.size())
            {
                throw std::length_error("BSSurfaceGeneral: wrong pole vector length.");
            }
            m_poles = std::move(poles);
        }

        /**
         * @brief Retrieves the boundary values of the surface in both U and V directions.
         * 
         * @return std::array<T, 4> The boundary values, where the first two elements are the boundaries in the U direction and the last two in the V direction.
         */
        virtual auto bounds() const -> std::array<T, 4> override
        {
            return m_bounds;
        }

        /**
         * @brief Changes the boundary values of the surface in the U direction.
         * 
         * This method updates the boundary values for the U direction and adjusts the associated knot vector accordingly.
         * 
         * @param k1 The new lower boundary value in the U direction.
         * @param k2 The new upper boundary value in the U direction.
         */
        auto changeUBounds(T k1, T k2) -> void
        {
            change_bounds(k1, k2, m_knotsFlatsU);
            m_bounds[0] = k1;
            m_bounds[1] = k2;
        }

        /**
         * @brief Changes the boundary values of the surface in the V direction.
         * 
         * Similar to changeUBounds, this method updates the boundary values for the V direction.
         * 
         * @param k1 The new lower boundary value in the V direction.
         * @param k2 The new upper boundary value in the V direction.
         */
        auto changeVBounds(T k1, T k2) -> void
        {
            change_bounds(k1, k2, m_knotsFlatsV);
            m_bounds[2] = k1;
            m_bounds[3] = k2;
        }

        /**
         * @brief Increases the degree of the surface in the U direction.
         * 
         * This method increases the polynomial degree of the surface in the U direction by one. It recalculates the control points and updates the knot vector.
         */
        auto increaseDegreeU() -> void
        {
            points_vector<T, dim+rational> poles_new;
            std::vector<T> ku_new;
            std::vector<T> ku_old{m_knotsFlatsU};
            auto nv = nPolesV();
            auto nu = nPolesU();
            for (size_t i{}; i < nv; i++)
            {

                points_vector<T, dim+rational> poles{std::next(m_poles.begin(), i * nu), std::next(m_poles.begin(), i * nu + nu)};
                if (i != 0)
                {
                    m_knotsFlatsU = {ku_old};
                }
                increase_degree(m_knotsFlatsU, poles, m_degU, 1);
                if (i == 0)
                {
                    ku_new = std::move(m_knotsFlatsU);
                }
                poles_new.insert(poles_new.end(), poles.begin(), poles.end());
            }
            m_poles = std::move(poles_new);
            m_knotsFlatsU = std::move(ku_new);
            m_degU++;
        }

        /**
         * @brief Increases the degree of the surface in the V direction.
         * 
         * This method increases the polynomial degree of the surface in the V direction by first swapping U and V directions, increasing the degree in U, and then swapping back.
         */
        auto increaseDegreeV() -> void
        {
            invertUV();
            increaseDegreeU();
            invertUV();
        }
 
        /**
         * @brief Swaps the U and V directions of the surface.
         * 
         * This method inverts the U and V directions of the surface, including swapping control points, knot vectors, and degrees. The boundaries are also updated accordingly.
         */
        auto invertUV() -> void
        {
            invert_uv_poles(m_poles,nPolesV());
            std::swap(m_knotsFlatsU,m_knotsFlatsV);
            std::swap(m_degU,m_degV);
            m_bounds = {m_knotsFlatsU.front(), m_knotsFlatsU.back(), m_knotsFlatsV.front(), m_knotsFlatsV.back()};
        }

        /**
         * @brief Reverses the direction of the surface in the U direction.
         * 
         * This method reverses the direction of the surface in the U direction by reversing the knot vector and the corresponding control points.
         */
        auto reverseU() -> void
        {
            auto k1 = m_knotsFlatsU.front();
            auto k2 = m_knotsFlatsU.back();
            std::reverse(m_knotsFlatsU.begin(), m_knotsFlatsU.end());
            std::transform(m_knotsFlatsU.begin(), m_knotsFlatsU.end(),
                           m_knotsFlatsU.begin(),
                           [&](const auto k_) {
                               return k1 + k2 - k_;
                           });
                        auto nv = nPolesV();
            auto nu = nPolesU();
            for (size_t i{}; i < nv; i++)
            {
                std::reverse(std::next(m_poles.begin(), i * nu), std::next(m_poles.begin(), i * nu + nu));
            }
        }

        /**
         * @brief Reverses the direction of the surface in the V direction.
         * 
         * This method reverses the direction of the surface in the V direction by first swapping U and V directions, reversing U, and then swapping back.
         */
        auto reverseV() -> void
        {
            invertUV();
            reverseU();
            invertUV();
        }

        /**
         * @brief Trims the surface in the U direction between two specified parameters.
         * 
         * This method adjusts the surface by trimming it between two parameter values in the U direction. 
         * It effectively narrows the surface to the specified U parameter range, updating the control points 
         * and knot vector in the U direction to reflect the new boundaries.
         * 
         * @param u1 The lower boundary of the new trimmed surface in the U direction.
         * @param u2 The upper boundary of the new trimmed surface in the U direction.
         */
        auto trimU(T u1, T u2) -> void
        {

            points_vector<T,dim+rational> new_poles;
            std::vector<T> new_knotsUFlats;
            for(size_t j{}; j < nPolesV(); j++)
            {
                new_knotsUFlats = std::vector<T>{m_knotsFlatsU};
                points_vector<T,dim+rational> poles_loc(nPolesU());
                std::copy(
                    std::next(m_poles.begin(), j * nPolesU()),
                    std::next(m_poles.begin(), (j+1) * nPolesU()),
                    poles_loc.begin()
                );
                trim(m_degU, new_knotsUFlats, poles_loc, u1, u2);
                new_poles.insert(new_poles.end(), poles_loc.begin(), poles_loc.end());
            }
            new_poles.shrink_to_fit();
            std::swap(m_knotsFlatsU, new_knotsUFlats);
            std::swap(m_poles, new_poles);

            m_bounds[0] = u1; 
            m_bounds[1] = u2;
        }

        /**
         * @brief Trims the surface in the V direction between two specified parameters.
         * 
         * This method adjusts the surface by trimming it between two parameter values in the V direction. 
         * It does this by first swapping the U and V directions, trimming in U (now representing V), 
         * and then swapping back. The control points and knot vector in the V direction are updated accordingly.
         * 
         * @param v1 The lower boundary of the new trimmed surface in the V direction.
         * @param v2 The upper boundary of the new trimmed surface in the V direction.
         */
        auto trimV(T v1, T v2) -> void
        {
            invertUV();
            trimU(v1, v2);
            invertUV();
        }

        /**
         * Creates an isoparametric curve on the surface at a specified parameter value in the U-direction.
         *
         * @param u The parameter value in the U-direction at which the isoparametric curve is to be created.
         * @return An isoparametric curve as a B-spline curve.
         * @throws OutOfBoundsSurfaceIsoU<T> if the parameter u is out of the defined bounds of the surface.
         * 
         */
        auto isoU(T u) const
        {
            auto [u1, u2, v1, v2] = this->bounds();
            if (u < u1 - knot_eps<T> || u > u2 + knot_eps<T>)
                throw OutOfBoundsSurfaceIsoU<T>(u,{u1, u2});
            // Proceed knot insertion on a copy, using forward declaration
            using srfType = typename std::conditional<rational, BSSurfaceRational<T,dim>, BSSurface<T,dim>>::type;
            srfType srf{*this};
            using crvType = typename std::conditional<rational, BSCurveRational<T,dim>, BSCurve<T,dim>>::type;
            auto j = srf.insertKnotU(u, srf.degreeU());
            return crvType{srf.polesV(j), srf.knotsFlatsV(), srf.degreeV()};
        }

        /**
         * Creates an isoparametric curve on the surface at a specified parameter value in the V-direction.
         *
         * @param v The parameter value in the V-direction at which the isoparametric curve is to be created.
         * @return An isoparametric curve as a B-spline curve.
         * @throws OutOfBoundsSurfaceIsoV<T> if the parameter v is out of the defined bounds of the surface.
         * 
         */
        auto isoV(T v) const
        {
            auto [u1, u2, v1, v2] = this->bounds();
            if (v < v1 - knot_eps<T> || v > v2 + knot_eps<T>)
                throw OutOfBoundsSurfaceIsoV<T>(v,{v1, v2});
            // Proceed knot insertion on a copy, using forward declaration
            using srfType = typename std::conditional<rational, BSSurfaceRational<T,dim>, BSSurface<T,dim>>::type;
            srfType srf{*this};
            using crvType = typename std::conditional<rational, BSCurveRational<T,dim>, BSCurve<T,dim>>::type;
            auto i = srf.insertKnotV(v, srf.degreeV());
            return crvType{srf.polesU(i), srf.knotsFlatsU(), srf.degreeU()};
        }
    };

    template <typename T, size_t dim>
    class BSSurface : public BSSurfaceGeneral<T, dim, false>
    {
    public:
        using BSSurfaceGeneral<T, dim, false>::BSSurfaceGeneral;
        BSSurface(const BSSurfaceGeneral<T, dim, false> &s)  : BSSurfaceGeneral<T, dim, false>{s} {}
        virtual auto value(T u, T v, size_t du = 0, size_t dv = 0) const -> std::array<T, dim> override
        {
            if (u < this->bounds()[0] - knot_eps<T> || u > this->bounds()[1] + knot_eps<T>)
                throw OutOfBoundsSurfaceUEval<T>(u,this->boundsU());
            if (v < this->bounds()[2] - knot_eps<T> || v > this->bounds()[3] + knot_eps<T>)
                throw OutOfBoundsSurfaceVEval<T>(v,this->boundsV());

            return eval_value_decasteljau(u, v, this->knotsFlatsU(), this->knotsFlatsV(), this->poles(), this->degreeU(), this->degreeV(), du, dv);
        }
    };

    template <typename T, size_t dim>
    class BSSurfaceRational : public BSSurfaceGeneral<T, dim, true>
    {
    public:
        using BSSurfaceGeneral<T, dim, true>::BSSurfaceGeneral;
        BSSurfaceRational(const BSSurfaceGeneral<T, dim, true> &s)  : BSSurfaceGeneral<T, dim, true>{s} {}
        BSSurfaceRational(const points_vector<T, dim> &poles,
                         const std::vector<T> &weights,
                         const std::vector<T> &knots_flatsU,
                         const std::vector<T> &knots_flatsV,
                         size_t degU,
                         size_t degV) : BSSurfaceGeneral<T, dim, true>{merge_weights(poles, weights), knots_flatsU, knots_flatsV, degU, degV} {}
        BSSurfaceRational(const std::vector<std::array<T, dim>> &poles,
                        const std::vector<T> &weights,
                        const std::vector<T> &knotsU,
                        const std::vector<T> &knotsV,
                        const std::vector<size_t> &multsU,
                        const std::vector<size_t> &multsV,
                        size_t degU,
                        size_t degV) : BSSurfaceGeneral<T, dim, true>{merge_weights(poles, weights), flat_knots(knotsU, multsU), flat_knots(knotsV, multsV), degU, degV} {}
        virtual auto value(T u, T v, size_t du = 0, size_t dv = 0) const -> std::array<T, dim> override
        {
            if (u < this->bounds()[0] - knot_eps<T> || u > this->bounds()[1] + knot_eps<T>)
                throw OutOfBoundsSurfaceUEval<T>(u,this->boundsU());
            if (v < this->bounds()[2] - knot_eps<T> || v > this->bounds()[3] + knot_eps<T>)
                throw OutOfBoundsSurfaceVEval<T>(v,this->boundsV());

            if (du == 0 && dv == 0)
            {
                return weight_projection(eval_value_decasteljau(u, v, this->knotsFlatsU(), this->knotsFlatsV(), this->poles(), this->degreeU(), this->degreeV(), du, dv));
            }
            else
            {
                throw std::runtime_error("not implemented");
            }
        }

        auto polesProjected() const -> points_vector<T, dim>
        {
            points_vector<T, dim> poles_(this->poles().size());
            std::transform(this->poles().begin(), this->poles().end(), poles_.begin(),
                           [](const auto &p)
                           { return weight_projection(p); });
            return poles_;
        }

        auto weights() const -> std::vector<T>
        {
            std::vector<T> weights_(this->poles().size());
            std::transform(this->poles().begin(), this->poles().end(), weights_.begin(),
                           [](const auto &p)
                           { return p.back(); });
            return weights_;
        }
    };

    template <typename T>
    class BSSfunction
    {
        BSSurface<T, 1> srf_;

    public:
        BSSfunction(const BSSfunction<T> &f) = default;
        BSSfunction(const BSSurface<T, 1> &crv) : srf_{crv} {}
        BSSfunction(
            const std::vector<T> &poles, 
            const std::vector<T> &knots_flatsU, 
            const std::vector<T> &knots_flatsV, 
            size_t degU, 
            size_t degV ) : srf_{BSCfunction<T>::poles_array(poles), knots_flatsU, knots_flatsV, degU, degV} {}
        BSSfunction(
            const std::vector<T> &poles, 
            const std::vector<T> &knotsU, 
            const std::vector<T> &knotsV,
            const std::vector<size_t> &multsU,
            const std::vector<size_t> &multsV,
            size_t degU, 
            size_t degV ) : srf_{BSCfunction<T>::poles_array(poles), knotsU, knotsV, multsU, multsV, degU, degV} {}
        auto operator()(T u, T v, size_t du = 0, size_t dv = 0) const -> T
        {
            return srf_.value(u, v, du, dv)[0];
        }
        auto bounds() -> std::array<T, 4> { return srf_.bounds(); }
        auto basisSurface() const -> const BSSurface<T,1> & {return srf_;}
    };

    using BSSurface2d_f = BSSurface<float, 2>;
    using BSSurface3d_f = BSSurface<float, 3>;
    using BSSurface2d_d = BSSurface<double, 2>;
    using BSSurface3d_d = BSSurface<double, 3>;
    using BSSurfaceRational2d_f = BSSurfaceRational<float, 2>;
    using BSSurfaceRational3d_f = BSSurfaceRational<float, 3>;
    using BSSurfaceRational2d_d = BSSurfaceRational<double, 2>;
    using BSSurfaceRational3d_d = BSSurfaceRational<double, 3>;

} // namespace gbs