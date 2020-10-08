#pragma once
#include <gbslib/bscanalysis.h>
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
    auto make_actor(const points_vector<T, dim> &pts,double pt_size=5.,bool render_as_sphere=true,double *col = {0.3,0.3,0.3}) -> vtkSmartPointer<vtkActor>
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

// class vtkShaderCallback : public vtkCommand
// {
// public:
//   static vtkShaderCallback *New()
//     { return new vtkShaderCallback; }
//   vtkRenderer *Renderer;
//   virtual void Execute(vtkObject *, unsigned long, void*cbo)
//     {
//     vtkOpenGLHelper *cellBO = reinterpret_cast<vtkOpenGLHelper*>(cbo);
//     cellBO->Program->SetUniformi("StipplePattern", 0xCC);
//     }

//   vtkShaderCallback() { this->Renderer = 0; }
// };



    template<typename T, size_t dim>
    auto make_actor(const points_vector<T,dim> &pts,const points_vector<T,dim> &poles) -> vtkSmartPointer<vtkAssembly>
    {
        auto colors = vtkSmartPointer<vtkNamedColors>::New();

        auto actor_crv = gbs::make_polyline(pts,colors->GetColor4d("Tomato").GetData());
        actor_crv->GetProperty()->SetLineWidth(3.f);

        auto ctrl_polygon = vtkSmartPointer<vtkAssembly>::New();

        vtkSmartPointer<vtkActor>  ctr_polygon_lines = gbs::make_polyline(poles,colors->GetColor4d("Black").GetData());
        auto ctr_polygon_dots = gbs::make_actor(poles,20.,true,colors->GetColor4d("Red").GetData()); 
        
        // std::ifstream in_geom("../render/stiple.geom");
        // std::string stiple((std::istreambuf_iterator<char>(in_geom)), 
        // std::istreambuf_iterator<char>());

        // std::ifstream in_dec("../render/stiple_dec.frag");
        // std::string stiple_dec((std::istreambuf_iterator<char>(in_dec)), 
        // std::istreambuf_iterator<char>());

        // std::ifstream in_impl("../render/stiple_impl.frag");
        // std::string stiple_impl((std::istreambuf_iterator<char>(in_impl)), 
        // std::istreambuf_iterator<char>());


        // ctr_polygon_lines->GetShaderProperty()->SetGeometryShaderCode(stiple.data());
        // ctr_polygon_lines->GetShaderProperty()->AddFragmentShaderReplacement("//VTK::TCoord::Dec",true,stiple_dec.data(),false);
        // ctr_polygon_lines->GetShaderProperty()->AddFragmentShaderReplacement("//VTK::TCoord::Impl",true,stiple_impl.data(),false);
        // // ctr_polygon_lines->GetShaderProperty()->GetFragmentCustomUniforms()->SetUniformi("StipplePattern",0xCC);
        // int size[2] = {800,600};
        // // ctr_polygon_lines->GetShaderProperty()->GetGeometryCustomUniforms()->SetUniform2i("ViewportSize",size);
        // ctr_polygon_lines->GetMapper()->AddObserver()



        ctrl_polygon->AddPart( ctr_polygon_lines );
        ctrl_polygon->AddPart( ctr_polygon_dots );

        ctr_polygon_lines->GetProperty()->SetLineWidth(3.f);
        gbs::StippledLine(ctr_polygon_lines,0xAAAA, 20);


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
        ctrl_polygon->AddPart(polesActor);

        auto ctr_polygon_lines = make_lattice_lines(poles,nPolesU,3.f,0.3,colors->GetColor3d("Black").GetData());

        ctrl_polygon->AddPart(ctr_polygon_lines);

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
        std::vector<std::array<T,dim+1>> p{bsc.poles()};

        points_vector<T, dim> poles;
        std::vector<T> weights;
        separate_weights(p, poles, weights);

        auto crv_actor = make_actor(pts, poles);
        auto crv_actor_parts = crv_actor->GetParts();
        auto ctrl_polygon = vtkAssembly::SafeDownCast(crv_actor_parts->GetItemAsObject(1));
        if (crv_actor)
        {
            auto sph_set = vtkAssembly::SafeDownCast(ctrl_polygon->GetParts()->GetItemAsObject(1));
            if (sph_set)
            {
                auto w = weights.begin();

                auto col = vtkPropCollection::New();
                sph_set->GetActors(col);

                auto it = col->NewIterator();
                for (it->InitTraversal(); !it->IsDoneWithTraversal(); it->GoToNextItem())
                {
                    auto actor_ = vtkActor::SafeDownCast(it->GetCurrentObject());
                    if (actor_)
                    {
                        actor_->SetScale(actor_->GetScale()[0]*(*w));
                        w++;
                    }
                }

                col->Delete();
            }
        }

        return crv_actor;
    }

        template<typename T, size_t dim>
    auto make_actor(const BSSurface<T,dim> &srf) -> vtkSmartPointer<vtkAssembly>
    {
        vtkIdType nu = 100 * srf.nPolesU();
        vtkIdType nv = 100 * srf.nPolesV();
        auto pts = gbs::discretize(srf,nu,nv); //TODO: improve discretization
        auto poles = srf.poles();
        std::vector<std::array<vtkIdType ,3> > pts_tri;

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
        
        // auto nu_iso = 10;
        // auto nv_iso = 10;
        // auto pts_iso = gbs::discretize(srf,nu_iso,nv_iso); //TODO: improve discretization
        // auto colors = vtkSmartPointer<vtkNamedColors>::New();
        // auto iso_lines = make_lattice_lines(pts_iso,nu_iso,1.f,1.,colors->GetColor3d("Black").GetData());
        // srf_actor->AddPart(iso_lines);
        
        return srf_actor;
    }

    auto make_actor(vtkProp3D* p){return p;}

    // template <typename container>
    // auto make_actor(const container &lst_) -> vtkSmartPointer<vtkAssembly>
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