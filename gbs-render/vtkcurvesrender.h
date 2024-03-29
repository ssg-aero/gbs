#pragma once
#include <gbs/bscanalysis.h>
#include <gbs/bssanalysis.h>
#include <gbs/bssinterp.h>
#include <vtkSmartPointer.h>
#include <vtkActor.h>
#include <vtkPoints.h>
#include <vtkAssembly.h>
#include <vtkProp3DCollection.h>
#include <vtkCollectionIterator.h>
#include <vtkNamedColors.h>
// #include <vtkRenderWindow.h>
// #include <vtkRenderWindowInteractor.h>
// #include <vtkRenderer.h>
#include <vtkPolyData.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkShaderProperty.h>
#include <vtkUniforms.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkUnstructuredGrid.h>
#include <vtkDataSetMapper.h>
#include <vtkPolyDataNormals.h>
#include <vtkAxesActor.h>
#include <vtkOrientationMarkerWidget.h>
#include <vtkLineSource.h>

namespace gbs
{
/**
 * @brief The default color used for points.
 *
 * This color is a static array of three doubles representing RGB values,
 * with default values of 0.3 for each component.
 */
    static std::array<double, 3> default_pnt_col{0.3, 0.3, 0.3};

/**
 * @brief Converts a color name to a color array.
 *
 * Given a color name, this function returns a color array as an std::array
 * of three doubles representing RGB values.
 *
 * @param color_name The name of the color to convert.
 *
 * @return The RGB color as an std::array<double, 3>.
 */
    auto array_from_col(const auto& color_name)
    {
        auto colors = vtkSmartPointer<vtkNamedColors>::New();
        auto col = colors->GetColor4d(color_name);
        return std::array<double, 3>{col.GetRed(), col.GetGreen(), col.GetBlue()};
    }
/**
 * @brief Convert general point to vtk point aka double[3]
 * 
 * @tparam T 
 * @tparam dim 
 * @param pt 
 * @return std::array<double,3> 
 */
    template<std::floating_point T,size_t dim>
    auto make_vtkPoint(const std::array<T,dim> &pt) -> std::array<double,3>
    {
        std::array<double,3> x = {0.,0.,0.};
        for(auto i = 0 ; i < fmin(dim,3);i++) x[i] = pt[i];
        return x;
    }
/**
 * @brief Create a VTK array of points from a generic container
 * 
 * @tparam Container 
 * @param pts 
 * @return vtkSmartPointer<vtkPoints> 
 */
    template <typename Container>
    auto make_vtkPoints(const Container &pts) -> vtkSmartPointer<vtkPoints>
    {
            vtkSmartPointer<vtkPoints> points =
        vtkSmartPointer<vtkPoints>::New();
        std::for_each(pts.begin(), pts.end(), [&](const auto pt_) { points->InsertNextPoint(make_vtkPoint(pt_).data()); });
        return points;
    }
/**
 * @brief Create a VTK array of points from a generic container
 * 
 * @tparam _FwdIt 
 * @param _First 
 * @param _Last 
 * @return vtkSmartPointer<vtkPoints> 
 */
    template <typename _FwdIt>
    auto make_vtkPoints(_FwdIt _First, _FwdIt _Last) -> vtkSmartPointer<vtkPoints>
    {
            vtkSmartPointer<vtkPoints> points =
        vtkSmartPointer<vtkPoints>::New();
        std::for_each(_First, _Last, [&](const auto &pt_) { points->InsertNextPoint(make_vtkPoint(pt_).data()); });
        return points;
    }
/**
 * @brief Create a vtkActor representing a polyline composed by the point array
 * 
 * @param pts 
 * @param color
 * @return vtkSmartPointer<vtkActor> 
 */
    GBS_EXPORT auto make_polyline_(vtkPoints *pts,std::array<double,3> &color) -> vtkSmartPointer<vtkActor>;
/**
 * @brief Change lines render as dashed
 * 
 * @param actor 
 * @param lineStipplePattern 
 * @param lineStippleRepeat 
 */
    GBS_EXPORT auto StippledLine(vtkSmartPointer<vtkActor> &actor,
                  int lineStipplePattern = 0xFFFF,
                  int lineStippleRepeat = 1) -> void;
/**
 * @brief Scale actors assembly
 * 
 * @param s 
 * @param a 
 */
    GBS_EXPORT auto scale_parts(double s,vtkAssembly *a) ->void;
/**
 * @brief Create a polyline actor from a container of points and color coordinates.
 *
 * @tparam Container A container type that stores the points.
 * @param pts A container of points.
 * @param coor An array of three doubles that represent the color coordinates.
 * @return A vtkSmartPointer<vtkActor> object representing the polyline.
 */
    template <typename Container>
    auto make_polyline(const Container &pts,std::array<double,3> &color) -> vtkSmartPointer<vtkActor>
    {
        return make_polyline_(make_vtkPoints(pts),color);
    }
/**
 * @brief Create a lattice of lines using the input points and other parameters.
 *
 * @tparam T The type of the input points.
 * @tparam dim The dimensionality of the input points.
 * @param points_msh The input points as a `points_vector`.
 * @param nU The number of points to extract along the u-direction.
 * @param lineWidth The width of the lines to be created.
 * @param opacity The opacity of the lines to be created.
 * @param col The color of the lines to be created as an array of three double values.
 * @return A smart pointer to a `vtkAssembly` object that contains the lines.
 */
    template <std::floating_point T, size_t dim>
    auto make_lattice_lines(const points_vector<T, dim> &points_msh, size_t nU,float lineWidth,double opacity,double *col) -> vtkSmartPointer<vtkAssembly>
    {
        auto lines = vtkSmartPointer<vtkAssembly>::New();
        auto nV = points_msh.size() / nU;
        for (auto iv = 0; iv < nU; iv++)
        {
            auto pt_iso = extract_V(iv, points_msh, nU);
            auto actor = make_polyline(pt_iso,col);
            actor->GetProperty()->SetLineWidth(lineWidth);
            actor->GetProperty()->SetOpacity(opacity);
            lines->AddPart(actor);
        }
        for (auto iu = 0; iu < nV; iu++)
        {
            auto poles_iso = extract_U(iu, points_msh, nU);
            auto actor = make_polyline(poles_iso, col);
            actor->GetProperty()->SetLineWidth(lineWidth);
            actor->GetProperty()->SetOpacity(opacity);
            lines->AddPart(actor);
        }
        return lines;
    }
/**
 * @brief Creates a VTK actor from a vector of points
 * 
 * @tparam T The type of the point coordinates
 * @tparam dim The dimensionality of the points
 * @tparam Container The type of the container holding the points
 * @param pts The container holding the points
 * @param pt_size The size of the points
 * @param render_as_sphere Whether to render the points as spheres
 * @param col The color of the points
 * @return vtkSmartPointer<vtkActor> The VTK actor representing the points
 */
    template <std::floating_point T, size_t dim>
    auto make_actor(const points_vector<T, dim> &pts,double pt_size=5.,bool render_as_sphere=true, std::array<double,3> &col=default_pnt_col) -> vtkSmartPointer<vtkActor>
    {
        auto Points = make_vtkPoints(pts);
        vtkSmartPointer<vtkPolyData> pointsPolydata =
            vtkSmartPointer<vtkPolyData>::New();

        pointsPolydata->SetPoints(Points);

        vtkSmartPointer<vtkVertexGlyphFilter> vertexFilter =
            vtkSmartPointer<vtkVertexGlyphFilter>::New();
        vertexFilter->SetInputData(pointsPolydata);
        vertexFilter->Update();

        vtkSmartPointer<vtkPolyData> polydata =
            vtkSmartPointer<vtkPolyData>::New();
        polydata->ShallowCopy(vertexFilter->GetOutput());

        vtkSmartPointer<vtkPolyDataMapper> pointMapper =
            vtkSmartPointer<vtkPolyDataMapper>::New();
        pointMapper->SetInputData(polydata);

        vtkSmartPointer<vtkActor> pointActor =
            vtkSmartPointer<vtkActor>::New();
        pointActor->SetMapper(pointMapper);
        pointActor->GetProperty()->SetPointSize(pt_size);
        pointActor->GetProperty()->SetRenderPointsAsSpheres(render_as_sphere);
        pointActor->GetProperty()->SetColor(col.data());

        return pointActor;
    }
/**
 * Create a control polygon VTK assembly.
 *
 * @tparam T            Data type for control points (e.g., double).
 * @tparam dim          Dimension of the control point space (e.g., 2 or 3).
 *
 * @param poles         Control points.
 * @param col_lines     Color for control polygon lines.
 * @param col_poles     Color for control polygon dots.
 *
 * @return              A VTK assembly containing the control polygon.
 */
    template <std::floating_point T, size_t dim>
    auto make_ctrl_polygon(const points_vector<T, dim> &poles, std::array<double,3> &col_lines, std::array<double,3> &col_poles) -> vtkSmartPointer<vtkAssembly>
    {
        auto ctrl_polygon = vtkSmartPointer<vtkAssembly>::New();


        vtkSmartPointer<vtkActor>  ctr_polygon_lines = make_polyline(poles,col_lines);
        auto ctr_polygon_dots = make_actor(poles,20.,true,col_poles); 
        
        ctrl_polygon->AddPart( ctr_polygon_lines );
        ctrl_polygon->AddPart( ctr_polygon_dots );

        // ctr_polygon_lines->SetPickable(false);

        ctr_polygon_lines->GetProperty()->SetLineWidth(3.f);
        ctr_polygon_lines->GetProperty()->SetOpacity(0.3);
        ctr_polygon_dots->GetProperty()->SetOpacity(0.3);
        // StippledLine(ctr_polygon_lines,0xAAAA, 20);
        return ctrl_polygon;
    }
/**
 * @brief Creates a curvature visualization for a 2D or 3D curve using VTK.
 * 
 * @tparam T The floating-point type used for the coordinates and computations.
 * @tparam dim The dimension of the curve (2 or 3).
 * @param crv The curve object for which the curvature visualization is to be created.
 * @param scale The scale factor for the curvature visualization.
 * @param np The number of points on the curve used to compute the curvature.
 * @param log_scale If true, the curvature values are scaled logarithmically.
 * @return A VTK Assembly containing the visualization of the curvature.
 */
    template <std::floating_point T, size_t dim>
    auto make_curvature(const Curve<T, dim> &crv, T scale = 0.1, size_t np = 30, bool log_scale = true) -> vtkSmartPointer<vtkAssembly>
    {
        auto curv_actor = vtkSmartPointer<vtkAssembly>::New();
        if(dim!=2 && dim!=3) return curv_actor;

        auto u = uniform_distrib_params<T, dim>(crv, np);
        points_vector<T, dim> pts { make_points<T, dim>(crv, u, 0) };
        points_vector<T, dim> pts_e{pts};

        std::transform(
            std::execution::par,
            pts.begin(),
            pts.end(),
            u.begin(),
            pts_e.begin(),
            [&](const auto &p, const auto &u_) {
                auto c = norm(crv(u_,2));
                if(log_scale) c = std::log10(1.+c);
                return p - normal_direction(crv,u_) *c * scale;
            });

        auto Tomato = array_from_col("Tomato");
        auto Lime = array_from_col("Lime");
        auto actor_cu = make_polyline(pts_e, Lime);
        actor_cu->GetProperty()->SetLineWidth(0.5f);


        curv_actor->AddPart(actor_cu);


        for (auto i = 0; i < np; i++)
        {
            vtkSmartPointer<vtkLineSource> lineSource = 
                vtkSmartPointer<vtkLineSource>::New();
            lineSource->SetPoint1(make_vtkPoint<T,dim>(pts[i]).data());
            lineSource->SetPoint2(make_vtkPoint<T,dim>(pts_e[i]).data());
            lineSource->Update();

            vtkSmartPointer<vtkPolyDataMapper> mapper = 
                vtkSmartPointer<vtkPolyDataMapper>::New();
            mapper->SetInputConnection(lineSource->GetOutputPort());
            vtkSmartPointer<vtkActor> actor =
                vtkSmartPointer<vtkActor>::New();
            actor->SetMapper(mapper);
            actor->GetProperty()->SetLineWidth(2);
            actor->GetProperty()->SetColor(Tomato.data());
            curv_actor->AddPart(actor);
        }
        curv_actor->AddPart(actor_cu);

        return curv_actor;
    }

/**
 * @brief Create a vtkAssembly object from a set of curve discretization points and its control points.
 * 
 * @tparam T Floating point type of the input points and control points.
 * @tparam dim Dimension of the input points and control points.
 * @param pts Input points to create the polyline actor.
 * @param poles Control points to create the control polygon actor.
 * @return vtkSmartPointer<vtkAssembly> Returns a vtkAssembly object which contains both the polyline and the control polygon actors.
 */
    template<std::floating_point T, size_t dim>
    auto make_actor(const points_vector<T,dim> &pts,const points_vector<T,dim> &poles) -> vtkSmartPointer<vtkAssembly>
    {
 
        double Tomato[3] = {255./255.,   99./255.,   71./255};
        double Black[3] = {0.,0.,0.} ;
        double Red[3] = {1.,0.,0.} ;
        auto actor_crv = make_polyline(pts,Tomato);
        actor_crv->GetProperty()->SetLineWidth(3.f);
        // actor_crv->SetPickable(false);

        auto ctrl_polygon = make_ctrl_polygon(poles,Black,Red);

        auto crv_actor = vtkSmartPointer<vtkAssembly>::New();
        crv_actor->AddPart(actor_crv);
        crv_actor->AddPart(ctrl_polygon);

        return crv_actor;
    }


    // template<std::floating_point T, size_t dim>
    // auto make_actor(const BSCurve<T,dim> &bsc, std::array<double,3> col = {255./255.,   99./255.,   71./255} ) //-> vtkSmartPointer<vtkAssembly>
    // {
    //     auto pts = discretize(bsc,30,0.01); 
    //     return make_polyline(pts,col.data());
    //     // auto pts = discretize(bsc,300); 
    //     // auto poles = bsc.poles();

    //     // return make_actor(pts,poles);

    // }

    // template <std::floating_point T, size_t dim>
    // auto make_actor(const BSCurveRational<T, dim> &bsc, std::array<double,3>  col = {255./255.,   99./255.,   71./255} ) //-> vtkSmartPointer<vtkAssembly>
    // {
    //     auto pts = discretize(bsc,30,0.01); 
    //     return make_polyline(pts,col.data());
    //     // auto pts = discretize(bsc,300); 
    //     // auto poles = bsc.polesProjected();

    //     // return make_actor(pts, poles);

    // }
/**
 * Create a VTK actor for a curve.
 *
 * @tparam T The type of the curve control points.
 * @tparam dim The dimension of the curve (2 or 3).
 * @param bsc The curve to create an actor for.
 * @param col The color of the curve.
 * @param np The number of points to use when discretizing the curve.
 * @param dev The maximum deviation between the curve and its discretization.
 * @return A VTK actor representing the curve.
 */
    template <std::floating_point T, size_t dim>
    auto make_actor(const Curve<T, dim> &bsc, std::array<double,3>  col = {255./255.,   99./255.,   71./255}, size_t np = 100, T dev = 0.01 ) //-> vtkSmartPointer<vtkAssembly>
    {
        auto pts = discretize<T, dim>(bsc,np,dev); 
        auto actor = make_polyline(pts,col);
        actor->GetProperty()->SetLineWidth(3.f);
        return actor;
    }


    template <std::floating_point T, size_t dim,bool rational>
    struct crv_dsp
    {
        const BSCurveGeneral<T,dim,rational> *c;
        std::array<T,3> col_crv   = {1.,0.,0.};
        bool poles_on = false;
        std::array<T,3> col_poles = {0.,1.,0.};
        std::array<T,3> col_ctrl  = {0.,0.,0.};
        float line_width = 3.f;
        bool show_curvature=false;
    };

    template <std::floating_point T,size_t dim,bool rational>
    auto make_actor(const crv_dsp<T,dim,rational> &cd) -> vtkSmartPointer<vtkAssembly>
    {

        auto pts = discretize(*cd.c, 30, 0.01);
        std::array<double,3> col_crv{cd.col_crv[0],cd.col_crv[1],cd.col_crv[2]}; // issue with const*
        std::array<double,3> col_ctrl{cd.col_ctrl[0],cd.col_ctrl[1],cd.col_ctrl[2]};
        std::array<double,3> col_poles{cd.col_poles[0],cd.col_poles[1],cd.col_poles[2]};

        auto actor_crv = make_polyline(pts,col_crv);
        actor_crv->GetProperty()->SetLineWidth(cd.line_width);
        
        auto crv_actor = vtkSmartPointer<vtkAssembly>::New();
        crv_actor->AddPart(actor_crv);
        if(cd.poles_on)
        {
            crv_actor->AddPart(make_ctrl_polygon(cd.c->poles(),col_ctrl,col_poles));
        }
        if(cd.show_curvature)
        {
            crv_actor->AddPart(make_curvature(*cd.c));
        }

        return crv_actor;
    }

} // namespace gbs