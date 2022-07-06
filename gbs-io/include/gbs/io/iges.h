#include <iges/api/dll_iges.h>
#include <iges/api/all_api_entities.h>
#include <gbs/curves>
#include <gbs/surfaces>
#include <gbs/bscapprox.h>

namespace gbs
{

    template <typename T, size_t d, bool rational>
    void add_geom(const BSCurveGeneral<T, d,rational> &crv, DLL_IGES &model,const std::string &name = "")
    {
        DLL_IGES_ENTITY_126 nc(model, true);
        auto [u1,u2] = crv.bounds();
        nc.SetNURBSData(
            crv.poles().size(),
            crv.degree() + 1,
            crv.knotsFlats().data(),
            const_cast<double *>(&crv.poles().data()[0][0]),
            rational,
            u1,
            u2
        );
        if(name.size()) nc.SetLabel(name.c_str());
    }

    template <typename T, size_t d, bool rational>
    void add_geom(const BSSurfaceGeneral<T, d,rational> &srf, DLL_IGES &model,const std::string &name = "")
    {
        DLL_IGES_ENTITY_128 nc(model, true);
        auto [u1,u2,v1,v2] = srf.bounds();
        nc.SetNURBSData(
            srf.nPolesU(),
            srf.nPolesV(),
            srf.degreeU() + 1,
            srf.degreeV() + 1,
            srf.knotsFlatsU().data(),
            srf.knotsFlatsV().data(),
            const_cast<double *>(&srf.poles().data()[0][0]),
            rational,false,false,
            u1,u2,v1,v2
        );
        if(name.size()) nc.SetLabel(name.c_str());
    }

    template <typename T, bool rational>
    void add_geom(const BSCurveGeneral<T,3,rational> &crv,const ax1<T,3> &ax,T v1, T v2, DLL_IGES &model,const std::string &name = "")
    {
        auto [u1, u2] = crv.bounds();
        DLL_IGES_ENTITY_120 rev( model, true );
        DLL_IGES_ENTITY_110 axis( model, true );
        DLL_IGES_ENTITY_126 nc(model, true);
        nc.SetNURBSData(
            crv.poles().size(),
            crv.degree() + 1,
            crv.knotsFlats().data(),
            const_cast<double *>(&crv.poles().data()[0][0]),
            rational,
            u1,
            u2
        );
        // axis
        axis.SetLineStart( ax[0][0],ax[0][1],ax[0][2] );
        axis.SetLineEnd( ax[0][0]+ax[1][0],ax[0][1]+ax[1][1],ax[0][2]+ax[1][2] );
        rev.SetAxis( axis );
        rev.SetGeneratrix( nc );
        rev.SetAngles( v1, v2 );
        if(name.size()) nc.SetLabel(name.c_str());
    }

    template <typename T>
    void add_geom(const SurfaceOfRevolution<T> &srf, DLL_IGES &model,const std::string &name = "")
    {
        const BSCurve<T,2> *p_bsc = dynamic_cast<const BSCurve<T,2>*>(srf.basisCurve().get());
        
        std::unique_ptr<BSCurve<T,2>> pu_bsc;
        if(!p_bsc) {
            pu_bsc = std::make_unique<BSCurve<T,2>>( approx(*srf.basisCurve(),0.01,5,KnotsCalcMode::CHORD_LENGTH, 1000) );
            p_bsc = pu_bsc.get();
        }


        auto poles2d = p_bsc->poles();
        gbs::points_vector<T,3> poles3d(poles2d.size());

        gbs::Matrix4<T> M = srf.transformation();
        std::transform(
            poles2d.begin(),
            poles2d.end(),
            poles3d.begin(),
            [&M](const auto &pt2d)
            {
                return gbs::transformed(add_dimension(pt2d),M);
            }
        );

        BSCurve<T,3> bsc3d{
            poles3d,
            p_bsc->knotsFlats(),
            p_bsc->degree()
        };

        auto [u1,u2,v1,v2] = srf.bounds();

        add_geom(bsc3d,srf.axis(),v1,v2,model,name);

    }

}