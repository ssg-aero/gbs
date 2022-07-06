#pragma once
#include <rapidjson/rapidjson.h>
#include <rapidjson/document.h>
#include <rapidjson/stringbuffer.h>
#include <rapidjson/writer.h>
#include <gbs/curves>
#include <gbs/surfaces>
#include <fstream>

namespace gbs
{

    // compatible with most iges types
    enum class entity_type {
        BSCfunction = 25,
        SurfaceOfRevolution = 120,
        BSCurve = 125,
        BSCurveRational = 126,
        BSSurface = 127,
        BSSurfaceRational = 128,
        // CurveOffset = 130,
        CurveTrimmed = 136,
        CurveComposite = 137,
        Curve2dOffset = 138,
        CurveOnSurface = 142
    };
    
    auto make_json(const auto &v_begin, const auto &v_end, auto &allocator) -> rapidjson::Value
    {
        rapidjson::Value v_val{rapidjson::kArrayType};
        std::for_each(
            v_begin, v_end,
            [&](auto value){v_val.PushBack(value,allocator);}
        );
        return v_val;
    }

    template< typename T>
    auto make_json(const std::shared_ptr<T> &p_v, auto &allocator) -> rapidjson::Value
    {
        return make_json(p_v.get(), allocator);
    }

    template< typename T, size_t dim>
    auto make_json(const std::vector<std::array<T,dim>> &v, auto &allocator) -> rapidjson::Value
    {
        rapidjson::Value v_val{rapidjson::kArrayType};
        std::for_each(
            v.begin(), v.end(),
            [&](const auto &value){
                v_val.PushBack(make_json(value.begin(), value.end(),allocator),allocator);
            }
        );
        return v_val;
    }

    template< typename T, size_t dim, size_t n>
    auto make_json(const std::array<std::array<T,dim>, n> &v, auto &allocator) -> rapidjson::Value
    {
        rapidjson::Value v_val{rapidjson::kArrayType};
        std::for_each(
            v.begin(), v.end(),
            [&](const auto &value){
                v_val.PushBack(make_json(value.begin(), value.end(),allocator),allocator);
            }
        );
        return v_val;
    }

    template< typename T, size_t dim>
    auto make_json(const BSCurve<T,dim> &crv, auto &allocator) -> rapidjson::Value
    {
        rapidjson::Value crv_val;
        crv_val.SetObject();
        rapidjson::Value deg_val{crv.degree()};
        rapidjson::Value dim_val{dim};
        rapidjson::Value type_val{static_cast<int>(entity_type::BSCurve)};
        auto [knots, mults] = knots_and_mults(crv.knotsFlats());
        auto knots_val = make_json(knots.begin(),knots.end(),allocator);
        auto mults_val = make_json(mults.begin(),mults.end(),allocator);
        auto poles_val = make_json(crv.poles(),allocator);

        crv_val.AddMember( "type"   ,type_val, allocator);
        crv_val.AddMember( "degree" ,deg_val, allocator);
        crv_val.AddMember( "dim"    ,dim_val, allocator);
        crv_val.AddMember( "knots"  ,knots_val, allocator);
        crv_val.AddMember( "mults"  ,mults_val, allocator);
        crv_val.AddMember( "poles"  ,poles_val, allocator);

        return crv_val;
    }

    template< typename T, size_t dim>
    auto make_json(const BSCurveRational<T,dim> &crv, auto &allocator) -> rapidjson::Value
    {
        rapidjson::Value crv_val;
        crv_val.SetObject();
        rapidjson::Value deg_val{crv.degree()};
        rapidjson::Value dim_val{dim};
        rapidjson::Value type_val{static_cast<int>(entity_type::BSCurveRational)};
        auto [knots, mults] = knots_and_mults(crv.knotsFlats());
        auto knots_val = make_json(knots.begin(),knots.end(),allocator);
        auto mults_val = make_json(mults.begin(),mults.end(),allocator);
        auto poles_val   = make_json(crv.polesProjected(),allocator);
        auto weights     = crv.weights();
        auto weights_val = make_json(weights.begin(),weights.end(),allocator);

        crv_val.AddMember( "type"    ,type_val, allocator);
        crv_val.AddMember( "degree"  ,deg_val, allocator);
        crv_val.AddMember( "dim"     ,dim_val, allocator);
        crv_val.AddMember( "knots"   ,knots_val, allocator);
        crv_val.AddMember( "mults"   ,mults_val, allocator);
        crv_val.AddMember( "poles"   ,poles_val, allocator);
        crv_val.AddMember( "weights" ,weights_val, allocator);

        return crv_val;
    }

    template< typename T>
    auto make_json(const BSCfunction<T> &f, auto &allocator) -> rapidjson::Value
    {
        rapidjson::Value f_val;
        f_val.SetObject();
        rapidjson::Value type_val{static_cast<int>(entity_type::BSCfunction)};
        auto crv_val = make_json<T,1>(f.basisCurve(), allocator);

        f_val.AddMember( "type"    ,type_val, allocator);
        f_val.AddMember( "curve"   ,crv_val, allocator);

        return f_val;
    }

    template< typename T>
    auto make_json(const gbs::CurveOffset<T, 2,gbs::BSCfunction<T>> & offset, auto &allocator) -> rapidjson::Value
    {
        rapidjson::Value offset_val;
        offset_val.SetObject();
        rapidjson::Value type_val{static_cast<int>(entity_type::Curve2dOffset)};
        rapidjson::Value dim_val{2};
        offset_val.AddMember( "dim"     ,dim_val, allocator);
        auto crv_val = make_json<T,2>(&offset.basisCurve(), allocator);
        auto f_val   = make_json(offset.offset(), allocator);

        offset_val.AddMember( "type"     ,type_val, allocator);
        offset_val.AddMember( "curve"    ,crv_val, allocator);
        offset_val.AddMember( "function" ,f_val, allocator);

        return offset_val;
    }

    template< typename T , size_t dim>
    auto make_json(const CurveComposite<T,dim> &crv, auto &allocator) -> rapidjson::Value
    {
        rapidjson::Value cc_val;
        cc_val.SetObject();
        rapidjson::Value type_val{static_cast<int>(entity_type::CurveComposite)};
        cc_val.AddMember( "type"     ,type_val, allocator);
        rapidjson::Value dim_val{dim};
        cc_val.AddMember( "dim"     ,dim_val, allocator);
        rapidjson::Value v_val{rapidjson::kArrayType};
        std::for_each(
            crv.curves().begin(), crv.curves().end(),
            [&](const auto &crv_){
                v_val.PushBack(
                    make_json(crv_, allocator),
                    allocator
                );
            }
        );
        cc_val.AddMember( "curves"    ,v_val, allocator);
        return cc_val;
    }

    template< typename T , size_t dim>
    auto make_json(const CurveTrimmed<T,dim> &crv, auto &allocator) -> rapidjson::Value
    {
        rapidjson::Value ct_val;
        ct_val.SetObject();
        rapidjson::Value type_val{static_cast<int>(entity_type::CurveTrimmed)};
        rapidjson::Value dim_val{dim};
        auto [u1, u2] = crv.bounds();
        rapidjson::Value u1_val{u1};
        rapidjson::Value u2_val{u2};
        ct_val.AddMember( 
            "type",
            type_val, 
            allocator
        );
        ct_val.AddMember( 
            "dim",
            dim_val,
            allocator
        );
        ct_val.AddMember( 
            "curve",
            make_json( crv.basisCurve(), allocator), 
            allocator
        );
        ct_val.AddMember(
            "u1",
            u1_val,
            allocator
        );
        ct_val.AddMember(
            "u2",
            u2_val,
            allocator
        );
        return ct_val;
    }


    template< typename T , size_t dim>
    auto make_json(const Curve<T,dim> *crv, auto &allocator) -> rapidjson::Value
    {
        if(dynamic_cast<const BSCurve<T,dim>*>(crv))
        {
            return make_json(*static_cast<const BSCurve<T,dim>*>(crv),allocator);
        }
        if(dynamic_cast<const BSCurveRational<T,dim>*>(crv))
        {
            return make_json(*static_cast<const BSCurveRational<T,dim>*>(crv),allocator);
        }
        if(dynamic_cast<const CurveOnSurface<T,dim>*>(crv))
        {
            return make_json(*static_cast<const CurveOnSurface<T,dim>*>(crv),allocator);
        }
        if(dynamic_cast<const CurveComposite<T,dim>*>(crv))
        {
            return make_json(*static_cast<const CurveComposite<T,dim>*>(crv),allocator);
        }
        if(dynamic_cast<const CurveTrimmed<T,dim>*>(crv))
        {
            return make_json(*static_cast<const CurveTrimmed<T,dim>*>(crv),allocator);
        }
        rapidjson::Value null_val;
        return null_val;
    }

    template< typename T>
    auto make_json(const Curve<T,2> *crv, auto &allocator) -> rapidjson::Value
    {
        if(dynamic_cast<const gbs::CurveOffset<T, 2,gbs::BSCfunction<T>>*>(crv))
        {
            return make_json(*static_cast<const gbs::CurveOffset<T, 2,gbs::BSCfunction<T>>*>(crv),allocator);
        }
        return make_json<T,2>(crv, allocator);
    }

    template <typename T >
    auto make_json(const SurfaceOfRevolution<T> &srf, auto &allocator) -> rapidjson::Value
    {
        rapidjson::Value srf_val;
        srf_val.SetObject();
        rapidjson::Value dim_val{2};
        rapidjson::Value type_val{static_cast<int>(entity_type::SurfaceOfRevolution)};
        auto ax2_val = make_json<T,3,3>(srf.axis2(), allocator);
        auto [u1, u2, th1, th2] = srf.bounds();
        rapidjson::Value th1_val{th1};
        rapidjson::Value th2_val{th2};
        auto crv_val = make_json(srf.basisCurve(), allocator);

        srf_val.AddMember( "type"    ,type_val, allocator);
        srf_val.AddMember( "dim"     ,dim_val, allocator);
        srf_val.AddMember("axis2",ax2_val, allocator);
        srf_val.AddMember("theta1",th1_val, allocator);
        srf_val.AddMember("theta2",th2_val, allocator);
        srf_val.AddMember("curve",crv_val, allocator);

        return srf_val;
    }

    template< typename T, size_t dim>
    auto make_json(const BSSurface<T,dim> &srf, auto &allocator) -> rapidjson::Value
    {
        rapidjson::Value srf_val;
        srf_val.SetObject();
        rapidjson::Value degU_val{srf.degreeU()};
        rapidjson::Value degV_val{srf.degreeV()};
        rapidjson::Value dim_val{dim};
        rapidjson::Value type_val{static_cast<int>(entity_type::BSSurface)};
        auto [knotsU, multsU] = knots_and_mults(srf.knotsFlatsU());
        auto [knotsV, multsV] = knots_and_mults(srf.knotsFlatsV());
        auto knotsU_val = make_json(knotsU.begin(),knotsU.end(),allocator);
        auto multsU_val = make_json(multsU.begin(),multsU.end(),allocator);
        auto knotsV_val = make_json(knotsV.begin(),knotsV.end(),allocator);
        auto multsV_val = make_json(multsV.begin(),multsV.end(),allocator);
        auto poles_val = make_json(srf.poles(),allocator);

        srf_val.AddMember( "type"    ,type_val, allocator);
        srf_val.AddMember( "degreeU" ,degU_val, allocator);
        srf_val.AddMember( "degreeV" ,degV_val, allocator);
        srf_val.AddMember( "dim"     ,dim_val, allocator);
        srf_val.AddMember( "knotsU"  ,knotsU_val, allocator);
        srf_val.AddMember( "multsU"  ,multsU_val, allocator);
        srf_val.AddMember( "knotsV"  ,knotsV_val, allocator);
        srf_val.AddMember( "multsV"  ,multsV_val, allocator);
        srf_val.AddMember( "poles"   ,poles_val, allocator);

        return srf_val;
    }

    template< typename T, size_t dim>
    auto make_json(const BSSurfaceRational<T,dim> &srf, auto &allocator) -> rapidjson::Value
    {
        rapidjson::Value srf_val;
        srf_val.SetObject();
        rapidjson::Value degU_val{srf.degreeU()};
        rapidjson::Value degV_val{srf.degreeV()};
        rapidjson::Value dim_val{dim};
        rapidjson::Value type_val{static_cast<int>(entity_type::BSSurfaceRational)};
        auto [knotsU, multsU] = knots_and_mults(srf.knotsFlatsU());
        auto [knotsV, multsV] = knots_and_mults(srf.knotsFlatsV());
        auto knotsU_val = make_json(knotsU.begin(),knotsU.end(),allocator);
        auto multsU_val = make_json(multsU.begin(),multsU.end(),allocator);
        auto knotsV_val = make_json(knotsV.begin(),knotsV.end(),allocator);
        auto multsV_val = make_json(multsV.begin(),multsV.end(),allocator);
        auto poles_val   = make_json(srf.polesProjected(),allocator);
        auto weights     = srf.weights();
        auto weights_val = make_json(weights.begin(),weights.end(),allocator);

        srf_val.AddMember( "type"     ,type_val, allocator);
        srf_val.AddMember( "degreeU"   ,degU_val, allocator);
        srf_val.AddMember( "degreeV"   ,degV_val, allocator);
        srf_val.AddMember( "dim"       ,dim_val, allocator);
        srf_val.AddMember( "knotsU"    ,knotsU_val, allocator);
        srf_val.AddMember( "multsU"    ,multsU_val, allocator);
        srf_val.AddMember( "knotsV"    ,knotsV_val, allocator);
        srf_val.AddMember( "multsV"    ,multsV_val, allocator);
        srf_val.AddMember( "poles"     ,poles_val, allocator);
        srf_val.AddMember( "weights"   ,weights_val, allocator);

        return srf_val;
    }

    template< typename T , size_t dim>
    auto make_json(const Surface<T,dim> *srf, auto &allocator) -> rapidjson::Value
    {
        if(dynamic_cast<const BSSurface<T,dim>*>(srf))
        {
            return make_json(*static_cast<const BSSurface<T,dim>*>(srf),allocator);
        }
        if(dynamic_cast<const BSSurfaceRational<T,dim>*>(srf))
        {
            return make_json(*static_cast<const BSSurfaceRational<T,dim>*>(srf),allocator);
        }
        rapidjson::Value null_val;
        return null_val;
    }

    template< typename T>
    auto make_json(const Surface<T,3> *srf, auto &allocator) -> rapidjson::Value
    {
        if(dynamic_cast<const SurfaceOfRevolution<T>*>(srf))
        {
            return make_json(*static_cast<const SurfaceOfRevolution<T>*>(srf),allocator);
        }
        return make_json<T,3>(srf, allocator);
    }

    template< typename T , size_t dim>
    auto make_json(const CurveOnSurface<T,dim> &crv, auto &allocator) -> rapidjson::Value
    {
        rapidjson::Value crv_val;
        crv_val.SetObject();
        rapidjson::Value dim_val{dim};
        rapidjson::Value type_val{static_cast<int>(entity_type::CurveOnSurface)};

        auto crv2d_val = make_json<T,2>(  &crv.basisCurve(),   allocator);
        auto srf_val   = make_json<T,dim>(&crv.basisSurface(), allocator);
        
        crv_val.AddMember( "type"     ,type_val,  allocator);
        crv_val.AddMember( "dim"      ,dim_val,   allocator);
        crv_val.AddMember( "curve2d"  ,crv2d_val, allocator);
        crv_val.AddMember( "surface"  ,srf_val,   allocator);

        return crv_val;
    }

    void write_js_doc(const rapidjson::Document &d, const char *fName)
    {
        rapidjson::StringBuffer buffer;
        rapidjson::Writer<rapidjson::StringBuffer> writer(buffer);
        d.Accept(writer);
        
        const char* output = buffer.GetString();

        std::ofstream myfile;
        myfile.open (fName);
        myfile << output;
        myfile.close();
    }
}