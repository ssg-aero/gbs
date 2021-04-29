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
        Line(const point<T,dim> &pnt, const point<T,dim> &dir) : ax_{pnt,dir/norm(dir)} {}
        Line(const ax1<T,dim> &ax) : Line{ax[0],ax[1]} {} // to force adim
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
                return {0.,0.};
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