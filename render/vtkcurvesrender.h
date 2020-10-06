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

namespace gbs
{
    template<typename T,size_t dim>
    auto make_vtkPoint(const std::array<T,dim> &pt) -> std::array<double,3>
    {
        std::array<double,3> x = {0.,0.,0.};
        for(auto i = 0 ; i < fmin(dim,3);i++) x[i] = pt[i];
        return x;
    }
    
    template <typename Container>
    auto make_vtkPoints(const Container &pts) -> vtkSmartPointer<vtkPoints>
    {
            vtkSmartPointer<vtkPoints> points =
        vtkSmartPointer<vtkPoints>::New();
        std::for_each(pts.begin(), pts.end(), [&](const auto pt_) { points->InsertNextPoint(make_vtkPoint(pt_).data()); });
        return points;
    }
    
    GBS_EXPORT auto make_polyline_(vtkPoints *pts,double *a) -> vtkSmartPointer<vtkActor>;

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

    template<typename T, size_t dim>
    auto make_actor(const points_vector<T,dim> &pts,const points_vector<T,dim> &poles) -> vtkSmartPointer<vtkAssembly>
    {
        auto colors = vtkSmartPointer<vtkNamedColors>::New();

        auto actor_crv = gbs::make_polyline(pts,colors->GetColor3d("Tomato").GetData());

        auto ctrl_polygon = vtkSmartPointer<vtkAssembly>::New();

        auto ctr_polygon_lines = gbs::make_polyline(poles,colors->GetColor3d("PaleGreen").GetData());
        auto ctr_polygon_dots = gbs::make_actor(poles,10.,true,colors->GetColor3d("Red").GetData()); 
        ctrl_polygon->AddPart( ctr_polygon_lines );
        ctrl_polygon->AddPart( ctr_polygon_dots );

        gbs::StippledLine(ctr_polygon_lines,0xAAAA, 2);


        auto crv_actor = vtkSmartPointer<vtkAssembly>::New();
        crv_actor->AddPart(actor_crv);
        crv_actor->AddPart(ctrl_polygon);

        return crv_actor;
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


    template <typename... Targs>
    auto plot(Targs... Fargs) -> void
    {
        vtkSmartPointer<vtkNamedColors> colors =
            vtkSmartPointer<vtkNamedColors>::New();

        // Setup render window, renderer, and interactor
        vtkSmartPointer<vtkRenderer> renderer =
            vtkSmartPointer<vtkRenderer>::New();
        vtkSmartPointer<vtkRenderWindow> renderWindow =
            vtkSmartPointer<vtkRenderWindow>::New();
        renderWindow->AddRenderer(renderer);
        vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
            vtkSmartPointer<vtkRenderWindowInteractor>::New();
        renderWindowInteractor->SetRenderWindow(renderWindow);

        auto tuple = std::tie(Fargs...);

        auto make_and_add_actor = [&](const auto &g){auto a = make_actor(g); renderer->AddActor(a);};

        tuple_for_each(tuple,make_and_add_actor);
                                        
        renderer->SetBackground(colors->GetColor3d("White").GetData());

        renderWindow->Render();
        renderWindowInteractor->Start();   
    }


} // namespace gbs