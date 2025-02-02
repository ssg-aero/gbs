#pragma once
#include <numbers>
#include <gbs/bscurve.h>
#include <gbs/transformpoints.h>
namespace gbs
{
    template <typename T, size_t dim>
    class Circle : public Curve<T, dim>
    {
        static_assert(dim == 2 || dim == 3, "Only 2d or 3d circle supported");
        ax2<T, dim> ax_;
        T r;

    public:
        Circle(const Circle<T,dim> &c) = default;
        Circle(T r, const ax2<T, dim> &ax) : r{r}, ax_{{ax[0], ax[1], r / norm(ax[2]) * ax[2] }} {} 
        Circle(T r, const ax1<T, dim> &ax) : r{r} {
            auto z = ax[1];
            auto i_min_coord_mag = std::distance(z.begin(), std::min_element(z.begin(), z.end()));
            point<T,dim> x{};
            x[i_min_coord_mag] = r;
            ax_ = { ax[0], ax[1], x };
        }

        virtual auto bounds() const -> std::array<T, 2> override
        {
            return {0, 2*std::numbers::pi_v<T>};
        }

        auto getAx() const noexcept -> const ax2<T, dim> &
        {
            return ax_;
        }

        auto getR() const noexcept -> T { return r; }
        /**
         * @brief Circle starting direction
         * 
         * @return const auto& 
         */
        const point<T,dim> & getX() const noexcept { return ax_[2];}
        /**
         * @brief Circle rotation axis
         * 
         * @return point<T,dim> 
         */
        const point<T,dim> & getZ() const noexcept { return ax_[1];}
        /**
         * @brief Circle center
         * 
         * @return point<T,dim> 
         */
        const point<T,dim> & getC() const noexcept { return ax_[0];}

    };

    template <typename T>
    class Circle2d : public Circle<T, 2>
    {
    public:
        Circle2d(const Circle2d<T> &c) = default;
        Circle2d(T r=1., const point<T, 2> &C={0,0}) : Circle<T, 2>{r, ax1<T, 2>{{{C}, point<T, 2>{0, 1}}}} { }; // {0,1} will set x to {1,0}
        virtual auto value(T u, size_t d = 0) const -> std::array<T, 2> override
        {

            u += d * 0.5 * std::numbers::pi_v<T>;

            point<T, 2> x = this->getX();
            rotate(x, u);
            if (d == 0)
            {
                x = x + this->getC();
            }
            return x;
        }
    };

    template <typename T>
    class Circle3d : public Circle<T, 3>
    {
    public:
        Circle3d(const Circle3d<T> &c) = default;
        Circle3d(T r, const ax2<T, 3> &ax) : Circle<T,3>{r, ax} { }
        Circle3d(T r=1, const ax1<T, 3> &ax={{0,0,0},{0,0,1}}) : Circle<T,3>{r, ax} { }
        virtual auto value(T u, size_t d = 0) const -> std::array<T, 3> override
        {

            u += d * 0.5 * std::numbers::pi_v<T>;

            point<T, 3> x = this->getX();
            rotate(x, u,this->getZ() );
            if (d == 0)
            {
                x = x + this->getC();
            }
            return x;
        }
    };

}