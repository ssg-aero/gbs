
#pragma once
#include <vector>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>
#include <vtkProperty.h>
#include <vtkAssembly.h>

#include "vtkcurvesrender.h"
#include "vtkSurfaceRender.h"

namespace gbs
{
    template <typename T>
    auto make_actor(const std::shared_ptr<T> &p_shr)
    {
        return make_actor(*p_shr);
    }

    template <typename T>
    auto make_actor(const std::vector<std::shared_ptr<T>> &lst_)
    {
        auto assembly_ = vtkSmartPointer<vtkAssembly>::New();
        std::ranges::for_each(lst_,
                             [&](const auto &c)
                             { assembly_->AddPart(make_actor(c)); });

        return assembly_;
    }

    inline auto make_actor(vtkProp3D *p) { return p; }

    template <typename T>
    auto make_actor(const std::vector<T> &lst_) -> vtkSmartPointer<vtkAssembly>
    {
        auto assembly_ = vtkSmartPointer<vtkAssembly>::New();
        std::ranges::for_each(lst_,
                             [&](const auto &c)
                             { assembly_->AddPart(make_actor(c)); });

        return assembly_;
    }

    /**
     * @brief : Add items to renderer and display a default VTK window
     *
     * @tparam Targs
     * @param Fargs
     */
    template <typename... Targs>
    auto plot(Targs... Fargs) -> vtkSmartPointer<vtkRenderWindow>
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
            vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New(); // like paraview

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

        auto make_and_add_actor = [&](const auto &g)
        {auto a = make_actor(g); renderer->AddActor(a); };

        tuple_for_each(tuple, make_and_add_actor);

        // renderer->AddActor(vtkSmartPointer<vtkAxesActor>::New());

        // renderer->SetBackground(colors->GetColor4d("White").GetData());
        renderer->SetBackground(0.9, 0.9, 0.95);

        renderer->ResetCamera();

        renderWindow->Render();
        renderWindowInteractor->Start();

        return renderWindow;
    }

} // namespace gbs