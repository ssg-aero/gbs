#pragma once
#include <iges/api/dll_iges.h>
#include <iges/api/all_api_entities.h>
#include <gbs/curves>
#include <gbs/surfaces>
#include <gbs/bscapprox.h>

#include <string>
#include <vector>

namespace gbs
{
    namespace iges_detail
    {
        /**
         * Writes a NURBS curve as an IGES entity 126. libIGES takes Cartesian
         * control points followed by the weight (x, y, z, w), always three
         * coordinates; gbs stores rational poles in homogeneous form
         * (w x, ..., w). 2D curves get z = 0. `scale` multiplies coordinates.
         */
        template <typename T, std::size_t dim, bool rational>
        void set_126(DLL_IGES_ENTITY_126 &e, const BSCurveGeneral<T, dim, rational> &c, T scale = T(1))
        {
            static_assert(dim >= 1 && dim <= 3, "IGES curves have at most 3 coordinates");
            std::vector<double> coeff;
            for (const auto &p : c.poles())
            {
                const T w = rational ? p[dim] : T(1);
                for (std::size_t k{}; k < dim; ++k)
                    coeff.push_back(double(p[k] / w * scale));
                for (std::size_t k{dim}; k < 3; ++k)
                    coeff.push_back(0.);
                if constexpr (rational)
                    coeff.push_back(double(w));
            }
            std::vector<double> knots(c.knotsFlats().begin(), c.knotsFlats().end());
            const auto [u1, u2] = c.bounds();
            e.SetNURBSData(int(c.poles().size()), int(c.order()), knots.data(), coeff.data(), rational, double(u1), double(u2));
        }

        /// Writes a NURBS surface as an IGES entity 128 (same pole convention as set_126; 2D surfaces get z = 0).
        template <typename T, std::size_t dim, bool rational>
        void set_128(DLL_IGES_ENTITY_128 &e, const BSSurfaceGeneral<T, dim, rational> &s, T scale = T(1))
        {
            static_assert(dim >= 1 && dim <= 3, "IGES surfaces have at most 3 coordinates");
            std::vector<double> coeff;
            for (const auto &p : s.poles())
            {
                const T w = rational ? p[dim] : T(1);
                for (std::size_t k{}; k < dim; ++k)
                    coeff.push_back(double(p[k] / w * scale));
                for (std::size_t k{dim}; k < 3; ++k)
                    coeff.push_back(0.);
                if constexpr (rational)
                    coeff.push_back(double(w));
            }
            std::vector<double> ku(s.knotsFlatsU().begin(), s.knotsFlatsU().end());
            std::vector<double> kv(s.knotsFlatsV().begin(), s.knotsFlatsV().end());
            const auto [u1, u2, v1, v2] = s.bounds();
            e.SetNURBSData(int(s.nPolesU()), int(s.nPolesV()), int(s.orderU()), int(s.orderV()), ku.data(), kv.data(),
                           coeff.data(), rational, false, false, double(u1), double(u2), double(v1), double(v2));
        }
    } // namespace iges_detail

    template <typename T, size_t d, bool rational>
    void add_geom(const BSCurveGeneral<T, d,rational> &crv, DLL_IGES &model, const std::string &name = "")
    {
        DLL_IGES_ENTITY_126 nc(model, true);
        iges_detail::set_126(nc, crv);
        if(name.size()) nc.SetLabel(name.c_str());
    }

    template <typename T, size_t d, bool rational>
    void add_geom(const BSSurfaceGeneral<T, d,rational> &srf, DLL_IGES &model, const std::string &name = "")
    {
        DLL_IGES_ENTITY_128 nc(model, true);
        iges_detail::set_128(nc, srf);
        if(name.size()) nc.SetLabel(name.c_str());
    }

    template <typename T, bool rational>
    void add_geom(const BSCurveGeneral<T,3,rational> &crv, const ax1<T,3> &ax, T v1, T v2, DLL_IGES &model, const std::string &name = "")
    {
        DLL_IGES_ENTITY_120 rev( model, true );
        DLL_IGES_ENTITY_110 axis( model, true );
        DLL_IGES_ENTITY_126 nc(model, true);
        iges_detail::set_126(nc, crv);
        // axis
        axis.SetLineStart( ax[0][0],ax[0][1],ax[0][2] );
        axis.SetLineEnd( ax[0][0]+ax[1][0],ax[0][1]+ax[1][1],ax[0][2]+ax[1][2] );
        rev.SetAxis( axis );
        rev.SetGeneratrix( nc );
        rev.SetAngles( v1, v2 );
        if(name.size()) nc.SetLabel(name.c_str());
    }

    template <typename T>
    void add_geom(const SurfaceOfRevolution<T> &srf, DLL_IGES &model, const std::string &name = "")
    {
        const BSCurve<T,2> *p_bsc = dynamic_cast<const BSCurve<T,2>*>(srf.basisCurve().get());
        
        std::unique_ptr<BSCurve<T,2>> pu_bsc;
        if(!p_bsc) {
            pu_bsc = std::make_unique<BSCurve<T,2>>( approx(*srf.basisCurve(),0.01,5,KnotsCalcMode::CHORD_LENGTH, 1000) );
            p_bsc = pu_bsc.get();
        }


        auto poles2d = p_bsc->poles();
        points_vector<T,3> poles3d(poles2d.size());

        Matrix4<T> M = srf.transformation();
        std::transform(
            poles2d.begin(),
            poles2d.end(),
            poles3d.begin(),
            [&M](const auto &pt2d)
            {
                return transformed(add_dimension(pt2d),M);
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

    template <typename T>
    class IgesWriter
    {
        DLL_IGES model_;
        public:
        IgesWriter() = default;
        void add_geometry(const BSCurve<T,2> &geom, const std::string &name = "")
        {
            add_geom<T,2>(geom, model_, name);
        }
        void add_geometry(const BSCurve<T,3> &geom, const std::string &name = "")
        {
            add_geom<T,3>(geom, model_, name);
        }
        void add_geometry(const BSCurveRational<T,2> &geom, const std::string &name = "")
        {
            add_geom<T,2>(geom, model_, name);
        }
        void add_geometry(const BSCurveRational<T,3> &geom, const std::string &name = "")
        {
            add_geom<T,3>(geom, model_, name);
        }
        void add_geometry(const BSSurface<T,2> &geom, const std::string &name = "")
        {
            add_geom<T,2>(geom, model_, name);
        }
        void add_geometry(const BSSurface<T,3> &geom, const std::string &name = "")
        {
            add_geom<T,3>(geom, model_, name);
        }
        void add_geometry(const BSSurfaceRational<T,2> &geom, const std::string &name = "")
        {
            add_geom<T,2>(geom, model_, name);
        }
        void add_geometry(const BSSurfaceRational<T,3> &geom, const std::string &name = "")
        {
            add_geom<T,3>(geom, model_, name);
        }
        void add_geometry(const SurfaceOfRevolution<T> &geom, const std::string &name = "")
        {
            add_geom<T>(geom, model_, name);
        }
        void write(const std::string &file_name, bool f_overwrite=true)
        {
            model_.Write(file_name.c_str(), f_overwrite);
        }
        /// Underlying libIGES model, for writers of other entities (gbs-io/iges_brep.h).
        DLL_IGES &model() noexcept { return model_; }
    };

}