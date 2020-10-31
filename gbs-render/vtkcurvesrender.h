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

    template<typename T, size_t dim>
    auto make_actor(const points_vector<T,dim> &pts,const points_vector<T,dim> &poles) -> vtkSmartPointer<vtkAssembly>
    {
        auto colors = vtkSmartPointer<vtkNamedColors>::New();

        auto actor_crv = gbs::make_polyline(pts,colors->GetColor4d("Tomato").GetData());
        actor_crv->GetProperty()->SetLineWidth(3.f);

        auto ctrl_polygon = vtkSmartPointer<vtkAssembly>::New();

        vtkSmartPointer<vtkActor>  ctr_polygon_lines = gbs::make_polyline(poles,colors->GetColor4d("Black").GetData());
        auto ctr_polygon_dots = gbs::make_actor(poles,20.,true,colors->GetColor4d("Red").GetData()); 
        
        ctrl_polygon->AddPart( ctr_polygon_lines );
        ctrl_polygon->AddPart( ctr_polygon_dots );

        ctr_polygon_lines->GetProperty()->SetLineWidth(3.f);
        ctr_polygon_lines->GetProperty()->SetOpacity(0.3);
        ctr_polygon_dots->GetProperty()->SetOpacity(0.3);
        // gbs::StippledLine(ctr_polygon_lines,0xAAAA, 20);


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

        return ugridActor;
    }

    template <typename T, size_t dim>
    auto make_actor(const points_vector<T, dim> &points_msh, const std::vector<std::array<vtkIdType, 3>> &pts_tri, const points_vector<T, dim> &poles,size_t nPolesU) -> vtkSmartPointer<vtkAssembly>
    {
        auto colors = vtkSmartPointer<vtkNamedColors>::New();

        auto srf_actor = vtkSmartPointer<vtkAssembly>::New();

        auto srf_msh_actor = make_actor(points_msh, pts_tri, colors->GetColor3d("Peacock").GetData());
        srf_actor->AddPart(srf_msh_actor);

        auto ctrl_polygon = vtkSmartPointer<vtkAssembly>::New();
        auto polesActor = gbs::make_actor(poles, 20., true, colors->GetColor4d("Red").GetData());
        polesActor->GetProperty()->SetOpacity(0.3);

        auto ctr_polygon_lines = make_lattice_lines(poles,nPolesU,3.f,0.3,colors->GetColor3d("Black").GetData());

        ctrl_polygon->AddPart(ctr_polygon_lines);
        ctrl_polygon->AddPart(polesActor);

        srf_actor->AddPart(ctrl_polygon);

        return srf_actor;
    }

    template<typename T, size_t dim>
    auto make_actor(const BSCurve<T,dim> &bsc) -> vtkSmartPointer<vtkAssembly>
    {
        auto pts = gbs::discretize(bsc,1000); //TODO: improve discretization
        auto poles = bsc.poles();

        return make_actor(pts,poles);

    }

    template <typename T, size_t dim>
    auto make_actor(const BSCurveRational<T, dim> &bsc) -> vtkSmartPointer<vtkAssembly>
    {
        auto pts = gbs::discretize(bsc, 1000); //TODO: improve discretization
        auto poles = bsc.polesProjected();

        return make_actor(pts, poles);

    }

    template <typename T, size_t dim>
    auto make_actor(const BSSurface<T, dim> &srf) -> vtkSmartPointer<vtkAssembly>
    {
        size_t n1 = 100 * srf.nPolesU();
        size_t n2 = 100 * srf.nPolesV();
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

        return  make_actor(pts, pts_tri, poles, srf.nPolesU());
    }

    template <typename T, size_t dim>
    auto make_actor(const BSSurfaceRational<T, dim> &srf) -> vtkSmartPointer<vtkAssembly>
    {
        size_t n1 = 100 * srf.nPolesU();
        size_t n2 = 100 * srf.nPolesV();
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

        auto srf_actor =  make_actor(pts,pts_tri,poles,srf.nPolesU());

        return srf_actor;
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


        renderWindowInteractor->SetInteractorStyle(style);

        auto tuple = std::tie(Fargs...);

        auto make_and_add_actor = [&](const auto &g){auto a = make_actor(g); renderer->AddActor(a);};

        tuple_for_each(tuple,make_and_add_actor);
                                        
        renderer->SetBackground(colors->GetColor4d("White").GetData());

        renderWindow->Render();
        renderWindowInteractor->Start();   

    }


} // namespace gbs