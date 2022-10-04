#pragma once
#include <limits>
#include <gbs/bscurve.h>
namespace gbs
{
    template <typename T, size_t dim>
    class Line : public Curve<T,dim>
    {
        ax1<T,dim> ax_;
        public:
        /**
         * @brief Construct a new Line objectbetwee 2 points
         * 
         * @param p1 
         * @param p2 
         */
        Line(const point<T,dim> &p1, const point<T,dim> &p2) : ax_{p1,(p2-p1)/norm(p2-p1)} {}
        /**
         * @brief Construct a new Line object from ax1
         * 
         * @param ax 
         */
        Line(const ax1<T,dim> &ax) : Line{ax[0],ax[0]+ax[1]} {} // to force adim

        Line(const Line<T,dim> &L) = default;

        virtual auto value(T u, size_t d = 0) const -> std::array<T, dim> override
        {
            switch (d)
            {
            case 0:
                return ax_[0] + u * ax_[1];
                break;
            case 1:
                return ax_[1];
                break;
            default:
                std::array<T, dim> res{};
                std::fill(res.begin(), res.end(),0);
                return res;
                break;
            }
        }
        virtual auto bounds() const -> std::array<T, 2> override
        {
            return {std::numeric_limits<T>::min(),std::numeric_limits<T>::max()};
        }

        auto getAx() const noexcept -> const ax1<T,dim>&
        {
            return ax_;
        }

        // virtual auto bounded() const -> bool {return false};

    };
    
}