#pragma once
#include "fromjson.h"

namespace gbs{
    template <typename T, size_t dim>
    auto bscurve_direct_mults(const rapidjson::Value &a) -> BSCurve<T, dim>
    {
        // std::cerr << "Curve name: " << a["name"].GetString() << std::endl;
        assert( a["type"].GetInt() == 125);
        assert( a["dim"].GetUint64() == dim);
        auto knots = make_vec<T>(a["knots"]);
        auto mults = make_vec<size_t>(a["mults"]);
        auto deg = get_val<size_t>(a["degree"]);
        auto poles = make_point_vec<T, dim>(a["poles"]);

        return BSCurve<T, dim>(poles, knots, mults, deg);
    }

}