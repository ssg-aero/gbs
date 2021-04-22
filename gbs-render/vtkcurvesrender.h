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
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
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
     * @brief Convert general point to vtk point aka double[3]
     * 
     * @tparam T 
     * @tparam dim 
     * @param pt 
     * @return std::array<double,3> 
     */
    template<typename T,size_t dim>
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
     * @param a 
     * @return vtkSmartPointer<vtkActor> 
     */
    GBS_EXPORT auto make_polyline_(vtkPoints *pts,double *a) -> vtkSmartPointer<vtkActor>;
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

    GBS_EXPORT auto scale_parts(double s,vtkAssembly *a) ->void;

    template <typename Container>
    auto make_polyline(const Container &pts,double *a) -> vtkSmartPointer<vtkActor>
    {
        return make_polyline_(make_vtkPoints(pts),a);
    }

    template <typename T, size_t dim>
    auto make_lattice_lines(const points_vector<T, dim> &points_msh, size_t nU,float lineWidth,double opacity,double *col) -> vtkSmartPointer<vtkAssembly>
    {
        auto lines = vtkSmartPointer<vtkAssembly>::New();
        auto nV = points_msh.size() / nU;
        for (auto iv = 0; iv < nU; iv++)
        {
            auto pt_iso = gbs::extract_V(iv, points_msh, nU);
            auto actor = gbs::make_polyline(pt_iso,col);
            actor->GetProperty()->SetLineWidth(lineWidth);
            actor->GetProperty()->SetOpacity(opacity);
            lines->AddPart(actor);
        }
        for (auto iu = 0; iu < nV; iu++)
        {
            auto poles_iso = gbs::extract_U(iu, points_msh, nU);
            auto actor = gbs::make_polyline(poles_iso, col);
            actor->GetProperty()->SetLineWidth(lineWidth);
            actor->GetProperty()->SetOpacity(opacity);
            lines->AddPart(actor);
        }
        return lines;
    }

    template <typename T, size_t dim>
    auto make_actor(const points_vector<T, dim> &pts,double pt_size=5.,bool render_as_sphere=true,double *col = std::array<double,3>{{0.3,0.3,0.3}}.data()) -> vtkSmartPointer<vtkActor>
    {
        auto Points = gbs::make_vtkPoints(pts);
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
        pointActor->GetProperty()->SetColor(col);

        return pointActor;
    }

    template <typename T, size_t dim>
    auto make_ctrl_polygon(const points_vector<T, dim> &poles,double *col_lines,double *col_poles) -> vtkSmartPointer<vtkAssembly>
    {
        auto ctrl_polygon = vtkSmartPointer<vtkAssembly>::New();


        vtkSmartPointer<vtkActor>  ctr_polygon_lines = gbs::make_polyline(poles,col_lines);
        auto ctr_polygon_dots = make_actor(poles,20.,true,col_poles); 
        
        ctrl_polygon->AddPart( ctr_polygon_lines );
        ctrl_polygon->AddPart( ctr_polygon_dots );

        // ctr_polygon_lines->SetPickable(false);

        ctr_polygon_lines->GetProperty()->SetLineWidth(3.f);
        ctr_polygon_lines->GetProperty()->SetOpacity(0.3);
        ctr_polygon_dots->GetProperty()->SetOpacity(0.3);
        // gbs::StippledLine(ctr_polygon_lines,0xAAAA, 20);
        return ctrl_polygon;
    }

    template <typename T, size_t dim>
    auto make_curvature(const Curve<T, dim> &crv, T scale = 1., size_t np = 30, bool log_scale = true) -> vtkSmartPointer<vtkAssembly>
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
                return p + normal_direction(crv,u_) *c * scale;
            });

        double Tomato[3] = {255./255.,   99./255.,   71./255};
        double Lime[3]   = {0.       ,         1.,        0.};
        auto actor_cu = gbs::make_polyline(pts_e, Lime);
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
            actor->GetProperty()->SetColor(Tomato);
            curv_actor->AddPart(actor);
        }
        curv_actor->AddPart(actor_cu);

        return curv_actor;
    }

    template<typename T, size_t dim>
    auto make_actor(const points_vector<T,dim> &pts,const points_vector<T,dim> &poles) -> vtkSmartPointer<vtkAssembly>
    {
 
        double Tomato[3] = {255./255.,   99./255.,   71./255};
        double Black[3] = {0.,0.,0.} ;
        double Red[3] = {1.,0.,0.} ;
        auto actor_crv = gbs::make_polyline(pts,Tomato);
        actor_crv->GetProperty()->SetLineWidth(3.f);
        // actor_crv->SetPickable(false);

        auto ctrl_polygon = make_ctrl_polygon(poles,Black,Red);

        auto crv_actor = vtkSmartPointer<vtkAssembly>::New();
        crv_actor->AddPart(actor_crv);
        crv_actor->AddPart(ctrl_polygon);

        return crv_actor;
    }

    template <typename T, size_t dim>
    auto make_actor(const points_vector<T, dim> &points_msh, const std::vector<std::array<vtkIdType, 3>> &pts_tri,double *col)
    {
        
        // auto pointActor = gbs::make_actor(pts,5.,true,colors->GetColor4d("Blue").GetData());

        vtkSmartPointer<vtkPolyData> ugrid =
            vtkSmartPointer<vtkPolyData>::New();
        ugrid->Allocate(pts_tri.size());

        std::for_each(pts_tri.begin(), pts_tri.end(),
                      [&ugrid](const auto &tri) { ugrid->InsertNextCell(VTK_TRIANGLE, 3, tri.data()); });

        ugrid->SetPoints( make_vtkPoints( points_msh ) );

        //TODO use surface's normals
        // Generate normals
        vtkSmartPointer<vtkPolyDataNormals> normalGenerator =
            vtkSmartPointer<vtkPolyDataNormals>::New();
        normalGenerator->SetInputData(ugrid);
        normalGenerator->ComputePointNormalsOn();
        normalGenerator->ComputeCellNormalsOff();
        normalGenerator->Update();

        /*
        // Optional settings
        normalGenerator->SetFeatureAngle(0.1);
        normalGenerator->SetSplitting(1);
        normalGenerator->SetConsistency(0);
        normalGenerator->SetAutoOrientNormals(0);
        normalGenerator->SetComputePointNormals(1);
        normalGenerator->SetComputeCellNormals(0);
        normalGenerator->SetFlipNormals(0);
        normalGenerator->SetNonManifoldTraversal(1);
        */

        ugrid = normalGenerator->GetOutput();

        vtkSmartPointer<vtkDataSetMapper> ugridMapper =
            vtkSmartPointer<vtkDataSetMapper>::New();
        ugridMapper->SetInputData(ugrid);

        vtkSmartPointer<vtkActor> ugridActor =
            vtkSmartPointer<vtkActor>::New();
        ugridActor->SetMapper(ugridMapper);
        ugridActor->GetProperty()->SetColor(col);
        // ugridActor->GetProperty()->EdgeVisibilityOn();
        ugridActor->GetProperty()->SetInterpolationToPhong();
        // ugridActor->GetProperty()->SetOpacity(0.7);

        return ugridActor;
    }

    template <typename T, size_t dim>
    auto make_actor(const points_vector<T, dim> &points_msh, const std::vector<std::array<vtkIdType, 3>> &pts_tri, const points_vector<T, dim> &poles,size_t nPolesU) -> vtkSmartPointer<vtkAssembly>
    {
        auto colors = vtkSmartPointer<vtkNamedColors>::New();

        auto srf_actor = vtkSmartPointer<vtkAssembly>::New();

        double Peacock[3] = { 51./255.,  161./255.,  201./255.};
        double Red[3] = { 1.,  0.,  0.};
        double Black[3] = { 0.,  0.,  0.};

        auto srf_msh_actor = make_actor(points_msh, pts_tri, Peacock);
        srf_actor->AddPart(srf_msh_actor);

        auto ctrl_polygon = vtkSmartPointer<vtkAssembly>::New();
        auto polesActor = gbs::make_actor(poles, 20., true, Red);
        polesActor->GetProperty()->SetOpacity(0.3);

        auto ctr_polygon_lines = make_lattice_lines(poles,nPolesU,3.f,0.3,Black);

        ctrl_polygon->AddPart(ctr_polygon_lines);
        ctrl_polygon->AddPart(polesActor);

        srf_actor->AddPart(ctrl_polygon);

        return srf_actor;
    }

    // template<typename T, size_t dim>
    // auto make_actor(const BSCurve<T,dim> &bsc, std::array<double,3> col = {255./255.,   99./255.,   71./255} ) //-> vtkSmartPointer<vtkAssembly>
    // {
    //     auto pts = gbs::discretize(bsc,30,0.01); 
    //     return make_polyline(pts,col.data());
    //     // auto pts = gbs::discretize(bsc,300); 
    //     // auto poles = bsc.poles();

    //     // return make_actor(pts,poles);

    // }

    // template <typename T, size_t dim>
    // auto make_actor(const BSCurveRational<T, dim> &bsc, std::array<double,3>  col = {255./255.,   99./255.,   71./255} ) //-> vtkSmartPointer<vtkAssembly>
    // {
    //     auto pts = gbs::discretize(bsc,30,0.01); 
    //     return make_polyline(pts,col.data());
    //     // auto pts = gbs::discretize(bsc,300); 
    //     // auto poles = bsc.polesProjected();

    //     // return make_actor(pts, poles);

    // }

    template <typename T, size_t dim>
    auto make_actor(const Curve<T, dim> &bsc, std::array<double,3>  col = {255./255.,   99./255.,   71./255}, size_t np = 100 ) //-> vtkSmartPointer<vtkAssembly>
    {
        auto pts = discretize<T, dim>(bsc,np,0.01); 
        auto actor = make_polyline(pts,col.data());
        actor->GetProperty()->SetLineWidth(3.f);
        return actor;
    }

    template <typename T, size_t dim>
    auto make_actor(const Surface<T, dim> &srf)// -> vtkSmartPointer<vtkAssembly>
    {
        size_t n1 = 200;
        size_t n2 = 200;
        auto pts = gbs::discretize(srf, n1, n2); //TODO: improve discretization

        std::vector<std::array<vtkIdType, 3>> pts_tri;

        vtkIdType nu = n1;
        vtkIdType nv = n2;

        std::array<vtkIdType, 3> tri;
        vtkIdType index;

        for (auto j = 0; j < nv - 1; j++)
        {
            for (auto i = 0; i < nu - 1; i++)
            {
                index = i + nu * j;
                pts_tri.push_back({index, index + 1, index + 1 + nu});
                pts_tri.push_back({index + 1 + nu, index + nu, index});
            }
        }

        double Peacock[3] = { 51./255.,  161./255.,  201./255.};
        return make_actor(pts, pts_tri, Peacock );
    }

    template <typename T, size_t dim>
    auto make_actor(const BSSurface<T, dim> &srf)// -> vtkSmartPointer<vtkAssembly>
    {
        size_t n1 = fmin(100 * srf.nPolesU(),1000);
        size_t n2 = fmin(100 * srf.nPolesV(),1000);
        auto pts = gbs::discretize(srf, n1, n2); //TODO: improve discretization
        auto poles = srf.poles();

        std::vector<std::array<vtkIdType, 3>> pts_tri;

        vtkIdType nu = n1;
        vtkIdType nv = n2;

        std::array<vtkIdType, 3> tri;
        vtkIdType index;

        for (auto j = 0; j < nv - 1; j++)
        {
            for (auto i = 0; i < nu - 1; i++)
            {
                index = i + nu * j;
                pts_tri.push_back({index, index + 1, index + 1 + nu});
                pts_tri.push_back({index + 1 + nu, index + nu, index});
            }
        }

        // return  make_actor(pts, pts_tri, poles, srf.nPolesU());
        double Peacock[3] = { 51./255.,  161./255.,  201./255.};
        return make_actor(pts, pts_tri, Peacock );
    }

    template <typename T, size_t dim>
    auto make_actor(const BSSurfaceRational<T, dim> &srf) //-> vtkSmartPointer<vtkAssembly>
    {
        size_t n1 = fmin(100 * srf.nPolesU(),1000);
        size_t n2 = fmin(100 * srf.nPolesV(),1000);
        auto pts = gbs::discretize(srf,n1,n2); //TODO: improve discretization
        auto poles = srf.polesProjected();
        
        std::vector<std::array<vtkIdType ,3> > pts_tri;

        vtkIdType nu = n1;
        vtkIdType nv = n2;

        std::array<vtkIdType ,3> tri;
        vtkIdType index;

        for (auto j = 0; j < nv - 1; j++)
        {
            for (auto i = 0; i < nu - 1; i++)
            {
                index = i + nu * j;
                pts_tri.push_back({index, index + 1, index + 1 + nu });
                pts_tri.push_back({index + 1 + nu, index + nu, index});
            }
        }

        // auto srf_actor =  make_actor(pts,pts_tri,poles,srf.nPolesU());
        // return srf_actor;

        auto colors = vtkSmartPointer<vtkNamedColors>::New();
        return make_actor(pts, pts_tri, colors->GetColor3d("Peacock").GetData());
    }

    template <typename T, size_t dim,bool rational>
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

    template <typename T,size_t dim,bool rational>
    auto make_actor(const crv_dsp<T,dim,rational> &cd) -> vtkSmartPointer<vtkAssembly>
    {

        auto pts = gbs::discretize(*cd.c, 30, 0.01);
        double col_crv[3] ={cd.col_crv[0],cd.col_crv[1],cd.col_crv[2]}; // issue with const*
        double col_ctrl[3] ={cd.col_ctrl[0],cd.col_ctrl[1],cd.col_ctrl[2]};
        double col_poles[3] ={cd.col_poles[0],cd.col_poles[1],cd.col_poles[2]};

        auto actor_crv = gbs::make_polyline(pts,col_crv);
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

    inline auto make_actor(vtkProp3D* p){return p;}

    template <typename T>
    auto make_actor(const std::vector<T> &lst_) -> vtkSmartPointer<vtkAssembly>
    {
        auto assembly_ = vtkSmartPointer<vtkAssembly>::New();
        std::for_each(lst_.begin(), lst_.end(),
                      [&](const auto &c) {
                                        assembly_->AddPart( gbs::make_actor(c) );});

        return assembly_;
    }

    /**
     * @brief : Add items to renderer and display a default VTK window
     * 
     * @tparam Targs 
     * @param Fargs 
     */
    template <typename... Targs>
    auto plot(Targs... Fargs) -> void
    {
        vtkSmartPointer<vtkNamedColors> colors =
            vtkSmartPointer<vtkNamedColors>::New();

        // Setup render window, renderer, and interactor
        vtkSmartPointer<vtkRenderer> renderer =
            vtkSmartPointer<vtkRenderer>::New();
        
        renderer->SetUseFXAA(true);

        vtkSmartPointer<vtkRenderWindow> renderWindow =
            vtkSmartPointer<vtkRenderWindow>::New();

        renderWindow->SetMultiSamples(0);

        renderWindow->AddRenderer(renderer);
        vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
            vtkSmartPointer<vtkRenderWindowInteractor>::New();
        renderWindowInteractor->SetRenderWindow(renderWindow);
        vtkSmartPointer<vtkInteractorStyleTrackballCamera> style =
            vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New(); //like paraview

        vtkSmartPointer<vtkAxesActor> axes =
            vtkSmartPointer<vtkAxesActor>::New();
        vtkSmartPointer<vtkOrientationMarkerWidget> widget =
            vtkSmartPointer<vtkOrientationMarkerWidget>::New();
        widget->SetOutlineColor(0.9300, 0.5700, 0.1300);
        widget->SetOrientationMarker(axes);
        widget->SetInteractor(renderWindowInteractor);
        widget->SetViewport(0.0, 0.0, 0.2, 0.2);
        widget->SetEnabled(1);
        widget->InteractiveOn();

        renderWindowInteractor->SetInteractorStyle(style);

        auto tuple = std::tie(Fargs...);

        auto make_and_add_actor = [&](const auto &g){auto a = make_actor(g); renderer->AddActor(a);};

        tuple_for_each(tuple,make_and_add_actor);

        // renderer->AddActor(vtkSmartPointer<vtkAxesActor>::New());
                                        
        // renderer->SetBackground(colors->GetColor4d("White").GetData());
        renderer->SetBackground(0.9,0.9,0.95);

        renderer->ResetCamera();

        renderWindow->Render();
        renderWindowInteractor->Start();   

    }


} // namespace gbs