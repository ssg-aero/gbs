#pragma once
#include <gbslib/bscurve.h>
#include <gbslib/transform.h>

namespace gbs
{
    template <typename T, size_t dim>
    auto build_ellipse(T radius1,T radius2,const std::array<T, dim> &center = std::array<T, dim>{} ) -> gbs::BSCurve<T, dim + 1>
    {
        std::vector<T> k = {0., 0., 0., 1. / 4., 1. / 4., 1. / 2., 1. / 2., 3. / 4., 3. / 4., 1., 1., 1.};
        auto wi = sqrt(2.) / 2.;

        std::vector<std::array<T, dim + 1>> poles(9);
        poles[0][0] = radius1;
        poles[0].back() = 1.;

        poles[1][0] = wi * radius1;
        poles[1][1] = wi * radius2;
        poles[1].back() = wi;

        poles[2][1] = radius2;
        poles[2].back() = 1;

        poles[3][0] = -wi * radius1;
        poles[3][1] = wi * radius2;
        poles[3].back() = wi;

        poles[4][0] = -radius1;
        poles[4].back() = 1;

        poles[5][0] = -wi * radius1;
        poles[5][1] = -wi * radius2;
        poles[5].back() = wi;

        poles[6][1] = -radius2;
        poles[6].back() = 1;

        poles[7][0] = wi * radius1;
        poles[7][1] = -wi * radius2;
        poles[7].back() = wi;

        poles[8][0] = radius1;
        poles[8].back() = 1;

        std::transform(
            poles.begin(), poles.end(),
            poles.begin(),
            [&center](auto p_) {
                translate(p_, center);
                return p_;
            });

        size_t p = 2;

        return gbs::BSCurve(poles, k, p);
    }

    template <typename T, size_t dim>
    auto build_circle(T radius,const std::array<T, dim> &center = std::array<T, dim>{} ) -> gbs::BSCurve<T, dim + 1>
    {
        return build_ellipse(radius,radius,center);
    }

    template <typename T, size_t dim>
    auto derivate(const BSCurve<T,dim> &crv) -> BSCurve<T,dim>
    {
        auto poles = crv.poles();
        auto knots = crv.knotsFlats();
        pointVector<T, dim> new_poles{crv.poles().size() - 1};
        auto p = crv.degree();
        auto n = new_poles.size();
        // std::vector<size_t> chunck(n);
        // std::generate(chunck.begin(), chunck.end(), [&, i = -1]() mutable {i++;return i; });
        // std::transform(
        //     std::execution::par,
        //     chunck.begin(), chunck.end(),
        //     new_poles.begin(),
        //     [&](const auto i) {
        //         return T(p) * (poles[i + 1] - poles[i]) / (knots[i + p + 1] - knots[i + 1]);
        //     });
        for(auto i = 0 ; i < n ; i++)
        {
            new_poles[i] = T(p) * (poles[i+1]-poles[i]) / (knots[i+p+1]-knots[i+1]);
        }
        knots.pop_back();
        knots.erase(knots.begin());
        return BSCurve<T,dim>(new_poles,knots,p-1);
    }

    template <typename T, size_t dim>
    auto integrate(const BSCurve<T,dim> &crv,std::array<T,dim> P0 = std::array<T,dim>{} ) -> BSCurve<T,dim>
    {
        auto poles = crv.poles();
        pointVector<T, dim> new_poles{poles.size() + 1};
        auto n = new_poles.size();
        auto knots = crv.knotsFlats();
        auto p = crv.degree()+1;
        knots.insert(knots.begin(),knots.front());
        knots.push_back(knots.back());
        new_poles[0] = P0;
        for(int i=1;i<n;i++)
        {
            new_poles[i] =  poles[i-1];
            new_poles[i] = new_poles[i]*(knots[p+i]-knots[i]);
            new_poles[i] = new_poles[i]/T(p);
            new_poles[i] = new_poles[i]+new_poles[i-1];
            // new_poles[i] =  poles[i-1] * (knots[p+i+1]-knots[i+1]) / T(p) + new_poles[i-1];
        }
        return BSCurve<T,dim>(new_poles,knots,p);
    }

} // namespace gbs