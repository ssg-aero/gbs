#include <gtest/gtest.h>
#include <gbs/bscbuild.h>
#include <gbs/bssbuild.h>
#include <gbs-io/tojson.h>
#include <gbs-io/fromjson.h>
#include <gbs-render/vtkGbsRender.h>
#include "rapidjson/stringbuffer.h"
#include <rapidjson/writer.h>

using namespace gbs;

#ifdef TEST_PLOT_ON
    const bool PLOT_ON = true;
#else
    const bool PLOT_ON = false;
#endif

TEST(tests_tojson, bsCurve)
{
    using T = double;
    const size_t dim =2;
    BSCurve<T,dim> crv{
        {
            {
                0.,0.
            },
            {
                0.5,0.
            },
            {
                0.5,1.
            },
            {
                1.,1.
            },
        },
        {
            0.,0.,0.,0.5,1.,1.,1.
        },
        2
    };

    rapidjson::Document document;
    document.SetObject();

    auto & crv_val = make_json(crv,document.GetAllocator());

    ASSERT_TRUE(crv_val["type"] == static_cast<int>(entity_type::BSCurve));
    ASSERT_TRUE(crv_val["degree"] == crv.degree());
    ASSERT_TRUE(crv_val["dim"] == dim);

    auto [knots, mults] = knots_and_mults(crv.knotsFlats());
    auto it_knots = knots.begin();
    for(const auto &k : crv_val["knots"].GetArray())
    {
        ASSERT_NEAR(k.GetDouble(),*it_knots, 1e-8);
        it_knots = std::next(it_knots);
    }
    auto it_mults = mults.begin();
    for(const auto &k : crv_val["mults"].GetArray())
    {
        ASSERT_NEAR(k.GetUint64(),*it_mults, 1e-8);
        it_mults = std::next(it_mults);
    }

    auto it_poles= crv.poles().begin();
    for(const auto &p : crv_val["poles"].GetArray())
    {
        auto pole = make_array<T,dim>(p);
        ASSERT_LT(distance(pole,*it_poles), 1e-8);
        it_poles = std::next(it_poles);
    }

    auto crv_ = make_curve<T,dim>(crv_val);
    ASSERT_LT(distance( crv_->begin(), crv.begin()), 1e-8);

    document.AddMember("crv",crv_val,document.GetAllocator());

    rapidjson::StringBuffer buffer;
    rapidjson::Writer<rapidjson::StringBuffer> writer(buffer);
    document.Accept(writer);
    
    const char* output = buffer.GetString();

    std::ofstream myfile;
    myfile.open ("../tests/out/tests_tojson_bscruve.json");
    myfile << output;
    myfile.close();
}


TEST(tests_tojson, bsCurveRational)
{
    using T = double;
    const size_t dim =2;
    auto crv = build_circle<T,2>(1.);

    rapidjson::Document document;
    document.SetObject();

    auto & crv_val = make_json(crv,document.GetAllocator());

    ASSERT_TRUE(crv_val["type"] == static_cast<int>(entity_type::BSCurveRational));
    ASSERT_TRUE(crv_val["degree"] == crv.degree());
    ASSERT_TRUE(crv_val["dim"] == dim);

    auto [knots, mults] = knots_and_mults(crv.knotsFlats());
    auto it_knots = knots.begin();
    for(const auto &k : crv_val["knots"].GetArray())
    {
        ASSERT_NEAR(k.GetDouble(),*it_knots, 1e-8);
        it_knots = std::next(it_knots);
    }
    auto it_mults = mults.begin();
    for(const auto &k : crv_val["mults"].GetArray())
    {
        ASSERT_NEAR(k.GetUint64(),*it_mults, 1e-8);
        it_mults = std::next(it_mults);
    }

    auto weights = crv.weights();
    auto it_weights = weights.begin();
    for(const auto &k : crv_val["weights"].GetArray())
    {
        ASSERT_NEAR(k.GetDouble(),*it_weights, 1e-8);
        it_weights = std::next(it_weights);
    }

    auto poles = crv.polesProjected();
    auto it_poles= poles.begin();
    for(const auto &p : crv_val["poles"].GetArray())
    {
        auto pole = make_array<T,dim>(p);
        ASSERT_LT(distance(pole,*it_poles), 1e-8);
        it_poles = std::next(it_poles);
    }

    document.AddMember("crv",crv_val,document.GetAllocator());

    rapidjson::StringBuffer buffer;
    rapidjson::Writer<rapidjson::StringBuffer> writer(buffer);
    document.Accept(writer);
    
    const char* output = buffer.GetString();

    std::ofstream myfile;
    myfile.open ("../tests/out/tests_tojson_bscruverational.json");
    myfile << output;
    myfile.close();

}

TEST(tests_tojson, bsSurface)
{
    using T = double;
    const size_t dim =3;
    BSSurface<T,dim> srf{
        {
            { 0., 0., 0. }, { 0.5, 0., 0. }, { 0.5, 1., 0. }, { 1., 1., 0. },
            { 0., 0., 1. }, { 0.5, 0., 1. }, { 0.5, 1., 1. }, { 1., 1., 1. },
        },
        {
            0.,0.,0.,0.5,1.,1.,1.
        },
        {
            0.,0.,1.,1.
        },
        2,
        1
    };

    rapidjson::Document document;
    document.SetObject();

    auto & srf_val = make_json(srf,document.GetAllocator());

    ASSERT_TRUE(srf_val["type"] == static_cast<int>(entity_type::BSSurface));
    ASSERT_TRUE(srf_val["degreeU"] == srf.degreeU());
    ASSERT_TRUE(srf_val["degreeV"] == srf.degreeV());
    ASSERT_TRUE(srf_val["dim"] == dim);

    auto [knotsU, multsU] = knots_and_mults(srf.knotsFlatsU());
    auto it_knotsU = knotsU.begin();
    for(const auto &k : srf_val["knotsU"].GetArray())
    {
        ASSERT_NEAR(k.GetDouble(),*it_knotsU, 1e-8);
        it_knotsU = std::next(it_knotsU);
    }
    auto it_multsU = multsU.begin();
    for(const auto &k : srf_val["multsU"].GetArray())
    {
        ASSERT_NEAR(k.GetUint64(),*it_multsU, 1e-8);
        it_multsU = std::next(it_multsU);
    }

    auto [knotsV, multsV] = knots_and_mults(srf.knotsFlatsV());
    auto it_knotsV = knotsV.begin();
    for(const auto &k : srf_val["knotsV"].GetArray())
    {
        ASSERT_NEAR(k.GetDouble(),*it_knotsV, 1e-8);
        it_knotsV = std::next(it_knotsV);
    }
    auto it_multsV = multsV.begin();
    for(const auto &k : srf_val["multsV"].GetArray())
    {
        ASSERT_NEAR(k.GetUint64(),*it_multsV, 1e-8);
        it_multsV = std::next(it_multsV);
    }

    auto it_poles= srf.poles().begin();
    for(const auto &p : srf_val["poles"].GetArray())
    {
        auto pole = make_array<T,dim>(p);
        ASSERT_LT(distance(pole,*it_poles), 1e-8);
        it_poles = std::next(it_poles);
    }

    document.AddMember("srf",srf_val,document.GetAllocator());

    rapidjson::StringBuffer buffer;
    rapidjson::Writer<rapidjson::StringBuffer> writer(buffer);
    document.Accept(writer);
    
    const char* output = buffer.GetString();

    std::ofstream myfile;
    myfile.open ("../tests/out/tests_tojson_bssurface.json");
    myfile << output;
    myfile.close();
}

TEST(tests_tojson, BSSurfaceRational)
{
    using T = double;
    const size_t dim =3;
    auto cir = build_circle<T,dim>(1.);
    auto ell = build_ellipse<T,dim>(1.,1.5,{0.,0.,1.});

    auto srf = loft(std::list<BSCurveGeneral<T,dim,true>*>{&cir, &ell});

    rapidjson::Document document;
    document.SetObject();

    auto & srf_val = make_json(srf,document.GetAllocator());

    ASSERT_TRUE(srf_val["type"] == static_cast<int>(entity_type::BSSurfaceRational));
    ASSERT_TRUE(srf_val["degreeU"] == srf.degreeU());
    ASSERT_TRUE(srf_val["degreeV"] == srf.degreeV());
    ASSERT_TRUE(srf_val["dim"] == dim);

    auto [knotsU, multsU] = knots_and_mults(srf.knotsFlatsU());
    auto it_knotsU = knotsU.begin();
    for(const auto &k : srf_val["knotsU"].GetArray())
    {
        ASSERT_NEAR(k.GetDouble(),*it_knotsU, 1e-8);
        it_knotsU = std::next(it_knotsU);
    }
    auto it_multsU = multsU.begin();
    for(const auto &k : srf_val["multsU"].GetArray())
    {
        ASSERT_NEAR(k.GetUint64(),*it_multsU, 1e-8);
        it_multsU = std::next(it_multsU);
    }

    auto [knotsV, multsV] = knots_and_mults(srf.knotsFlatsV());
    auto it_knotsV = knotsV.begin();
    for(const auto &k : srf_val["knotsV"].GetArray())
    {
        ASSERT_NEAR(k.GetDouble(),*it_knotsV, 1e-8);
        it_knotsV = std::next(it_knotsV);
    }
    auto it_multsV = multsV.begin();
    for(const auto &k : srf_val["multsV"].GetArray())
    {
        ASSERT_NEAR(k.GetUint64(),*it_multsV, 1e-8);
        it_multsV = std::next(it_multsV);
    }

    auto poles = srf.polesProjected();
    auto it_poles= poles.begin();
    for(const auto &p : srf_val["poles"].GetArray())
    {
        auto pole = make_array<T,dim>(p);
        ASSERT_LT(distance(pole,*it_poles), 1e-8);
        it_poles = std::next(it_poles);
    }

    auto weights = srf.weights();
    auto it_weights = weights.begin();
    for(const auto &k : srf_val["weights"].GetArray())
    {
        ASSERT_NEAR(k.GetDouble(),*it_weights, 1e-8);
        it_weights = std::next(it_weights);
    }

    document.AddMember("srf",srf_val,document.GetAllocator());

    rapidjson::StringBuffer buffer;
    rapidjson::Writer<rapidjson::StringBuffer> writer(buffer);
    document.Accept(writer);
    
    const char* output = buffer.GetString();

    std::ofstream myfile;
    myfile.open ("../tests/out/tests_tojson_bssurfacerational.json");
    myfile << output;
    myfile.close();

    if(PLOT_ON) plot(srf);  
}

TEST(tests_tojson, BSSurfaceOfRevolution)
{
    using T = double;
    auto cir = std::make_shared<BSCurveRational<T,2>>(build_circle<T,2>(0.1, {0., 1.}));

    ax2<T,3> axz{{
        {0.,0.,0.},
        {0.,0.,1.},
        {1.,0.,0.}
    }};
    SurfaceOfRevolution<T> srf(cir, axz);

    rapidjson::Document document;
    document.SetObject(); 

    auto & srf_val = make_json(srf,document.GetAllocator());

    document.AddMember("srf",srf_val,document.GetAllocator());

    rapidjson::StringBuffer buffer;
    rapidjson::Writer<rapidjson::StringBuffer> writer(buffer);
    document.Accept(writer);
    
    const char* output = buffer.GetString();

    std::ofstream myfile;
    myfile.open ("../tests/out/tests_tojson_surfacerevolution.json");
    myfile << output;
    myfile.close();

    if(PLOT_ON) plot(srf);  
}

TEST(tests_tojson, CurveOnSurface)
{
    using T = double;
    const size_t dim =3;
    auto cir = build_circle<T,dim>(1.);
    auto ell = build_ellipse<T,dim>(1.,1.5,{0.,0.,1.});

    auto srf = std::make_shared<BSSurfaceRational<T,dim>>(
        loft(std::list<BSCurveGeneral<T,dim,true>*>{&cir, &ell})
    );

    auto crv2d = std::make_shared<BSCurve<T,2>>(
        BSCurve<T,2>{
            {
                {0.,0.},{1.,1.}
            },
            {0., 0.,1.,1.},
            1
        }
    );

    CurveOnSurface<T,dim> crv{
        crv2d, srf
    };

    rapidjson::Document document;
    document.SetObject(); 

    auto & crv_val = make_json(crv,document.GetAllocator());

    document.AddMember("crv",crv_val,document.GetAllocator());

    rapidjson::StringBuffer buffer;
    rapidjson::Writer<rapidjson::StringBuffer> writer(buffer);
    document.Accept(writer);
    
    const char* output = buffer.GetString();

    std::ofstream myfile;
    myfile.open ("../tests/out/tests_tojson_curveonsurface.json");
    myfile << output;
    myfile.close();
    if(PLOT_ON)
        plot(srf,crv);
}

TEST(tests_tojson, Curve2dOffset)
{
    auto circle = gbs::build_circle<double, 2>(1.);
    auto p_circle = std::make_shared<gbs::BSCurveRational<double, 2>>(circle);
    auto f_offset = gbs::BSCfunction<double>(gbs::build_segment<double, 1>({-1.}, {-1.},true));
    gbs::CurveOffset2D<double, gbs::BSCfunction<double>> circle2{
        p_circle,
        f_offset};


    rapidjson::Document document;
    document.SetObject(); 

    auto & crv_val = make_json(circle2,document.GetAllocator());

    document.AddMember("crv",crv_val,document.GetAllocator());

    rapidjson::StringBuffer buffer;
    rapidjson::Writer<rapidjson::StringBuffer> writer(buffer);
    document.Accept(writer);
    
    const char* output = buffer.GetString();

    std::ofstream myfile;
    myfile.open ("../tests/out/tests_tojson_curve2doffset.json");
    myfile << output;
    myfile.close();

    
}