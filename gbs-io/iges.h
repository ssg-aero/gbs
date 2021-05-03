#include <iges/api/dll_iges.h>
#include <iges/api/all_api_entities.h>
#include <gbs/curves>
#include <gbs/surfaces>

namespace gbs
{

    template <typename T, size_t d, bool rational>
    void add_geom(const BSCurveGeneral<T, d,rational> &crv, DLL_IGES &model)
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
    }

    template <typename T, size_t d, bool rational>
    void add_geom(const BSSurfaceGeneral<T, d,rational> &srf, DLL_IGES &model)
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
    }

    template <typename T, bool rational>
    void add_geom(const BSCurveGeneral<T,3,rational> &crv,const ax1<T,3> &ax,T v1, T v2, DLL_IGES &model)
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
    }

    template <typename T>
    void add_geom(const SurfaceOfRevolution<T> &srf, DLL_IGES &model)
    {
        const BSCurve<T,2> *bsc = static_cast<const BSCurve<T,2>*>(srf.basisCurve().get());
        
        if(!bsc) return;

        auto poles2d = bsc->poles();
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
            bsc->knotsFlats(),
            bsc->degree()
        };

        auto [u1,u2,v1,v2] = srf.bounds();

        add_geom(bsc3d,srf.axis(),v1,v2,model);

        // DLL_IGES_ENTITY_120 rev( model, true );
        // DLL_IGES_ENTITY_110 axis( model, true );
        // DLL_IGES_ENTITY_126 nc(model, true);
        // nc.SetNURBSData(
        //     bsc3d.poles().size(),
        //     bsc3d.degree() + 1,
        //     bsc3d.knotsFlats().data(),
        //     const_cast<double *>(&bsc3d.poles().data()[0][0]),
        //     false,
        //     u1,
        //     u2
        // );
        // // axis
        // auto ax_ = srf.axis();
        // axis.SetLineStart( ax_[0][0],ax_[0][1],ax_[0][2] );
        // axis.SetLineEnd( ax_[0][0]+ax_[1][0],ax_[0][1]+ax_[1][1],ax_[0][2]+ax_[1][2] );
        // rev.SetAxis( axis );
        // rev.SetGeneratrix( nc );
        // rev.SetAngles( v1, v2 );

        // add_geom(bsc3d,model);
    }

}