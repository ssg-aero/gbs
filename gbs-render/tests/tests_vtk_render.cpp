#include <gtest/gtest.h>

#include <gbs/bscanalysis.h>
#include <gbs/bssanalysis.h>
#include <gbs/bscbuild.h>
#include <gbs/bscurve.h>
#include <gbs/bscapprox.h>

#include <gbs-occt/export.h>
#include <gbs-occt/curvesbuild.h>

#include <gbs-render/vtkcurvesrender.h>


using gbs::operator+;
using gbs::operator/;


TEST(tests_vtk_render, BSC)
{

    auto r1 = 3.;
    auto r2 = 1.;
    auto c = 25.;
    auto crv1 = gbs::build_ellipse<double, 3>(r1, r2, {r1, 0., 0.});
    crv1.trim(0.25, 0.5);
    crv1.reverse();

    crv1.changeBounds(0,r1);
   
    std::vector<std::array<double, 4>> poles2 = {{r1,r2,0.,1.},{0.5*(r1+c),r2,0.,1.},{c,r2,0.,1.}};
    std::vector<double>                    k2 = {r1,r1,r1,c,c,c};

    gbs::BSCurveRational3d_d crv2(poles2,k2,2);
    
    auto k1 = crv1.knotsFlats();
    auto poles1 = crv1.poles();
    k1.pop_back();
    k1.pop_back();
    k1.pop_back();
    k2.erase(k2.begin());
    poles1.pop_back();

    k1.insert(k1.end(), k2.begin(), k2.end());
    poles1.insert(poles1.end(), poles2.begin(), poles2.end());

    gbs::BSCurveRational3d_d crv3(poles1,k1,2); 
    crv3.changeBounds(0.,1.);

    std::vector<gbs::BSCurveRational3d_d> c_lsr{crv3};
    auto poles_bis = crv3.poles();
    poles_bis.back()[1]+=r2;

    gbs::points_vector_3d_d poles_non_rational(crv3.poles().size());
    std::transform(crv3.poles().begin(),crv3.poles().end(),poles_non_rational.begin(),[](const auto &p_){return gbs::weight_projection(p_);});
    
    gbs::plot(
        c_lsr,
        gbs::BSCurveRational3d_d(poles_bis,crv3.knotsFlats(),crv3.degree()),
    	gbs::BSCurve3d_d(poles_non_rational,crv3.knotsFlats(),crv3.degree())
        );

}

TEST(tests_vtk_render, BSC_option)
{
    std::vector<double> k1 = {0., 0., 0., 0., 1., 1., 1., 1.};
    std::vector<double> k2 = {0., 0., 0., 0., 1., 1., 1., 1.};
    std::vector<std::array<double,3> > poles1 =
    {
        {0.,1.,0.},
        {1.,2.,0.},
        {2.,2.,0.},
        {3.,0.,0.},
    };

    size_t p1 = 3;

    gbs::BSCurve3d_d c1(poles1,k1,p1);

    gbs::plot(
        gbs::crv_dsp<double,3,false>{
            .c =&c1,
            .col_crv = {1.,0.,0.},
            .poles_on = true,
            .col_poles = {0.,1.,0.},
            .col_ctrl = {0.,0.,0.},
            .show_curvature=true,
            } // c++20
        );

}

TEST(tests_vtk_render, dev)
{

    auto r1 = 3.;
    auto r2 = 1.;
    auto c = 25.;
    auto crv1 = gbs::build_ellipse<double, 3>(r1, r2, {r1, 0., 0.});
    crv1.trim(0.25, 0.5);
    crv1.reverse();

    crv1.changeBounds(0,r1);
   
    std::vector<std::array<double, 4>> poles2 = {{r1,r2,0.,1.},{0.5*(r1+c),r2,0.,1.},{c,r2,0.,1.}};
    std::vector<double>                    k2 = {r1,r1,r1,c,c,c};

    gbs::BSCurveRational3d_d crv2(poles2,k2,2);
    
    auto k1 = crv1.knotsFlats();
    auto poles1 = crv1.poles();
    k1.pop_back();
    k1.pop_back();
    k1.pop_back();
    k2.erase(k2.begin());
    poles1.pop_back();

    k1.insert(k1.end(), k2.begin(), k2.end());
    poles1.insert(poles1.end(), poles2.begin(), poles2.end());

    gbs::BSCurveRational3d_d crv3(poles1,k1,2); // thickenss law
    crv3.changeBounds(0.,1.);

    auto b1 = 60 / 180. * M_PI;
    auto b2 = -10 / 180. * M_PI;
    auto t1 = 1.;
    auto t2 = 1.;
    auto g  = 30 / 180. * M_PI;

    std::vector<std::array<double, 4>> poles4 = {{0.,0.,0.,1.},{0.33*t1*c*cos(b1),0.33*t1*c*sin(b1),0.,1.},{c*cos(g)-0.33*t2*c*cos(b2),c*sin(g)-0.33*t2*c*sin(b2),0.,1.},{c*cos(g),c*sin(g),0.,1.}};
    std::vector<double>                    k4 = {0.,0.,0.,0.5,1.,1.,1.};
    gbs::BSCurveRational3d_d crv4(poles4,k4,2); // mean line

    //curves have the same start/stop
    // attention dansle cas général cela ne marche pas

    // curves uniformization
    auto k3 = crv3.knotsFlats();
    std::for_each(std::next(k4.begin(),3),std::next(k4.end(),-3),[&](const auto k_){crv3.insertKnot(k_);});
    std::for_each(std::next(k3.begin(),3),std::next(k3.end(),-3),[&](const auto k_){crv4.insertKnot(k_);});

    std::vector<double> u(crv3.poles().size());

    u[0] = 0.;
    u[1] = 0.;
    std::transform(std::next(crv3.poles().begin(), 2), crv3.poles().end(), std::next(u.begin(), 2),
                   [&](const auto &pt_) {
                       return gbs::extrema_curve_point(crv3, gbs::weight_projection(pt_), 1e-6)[0];
                   });

    std::vector<std::array<double, 4>> poles5(crv3.poles().size());
    // std::transform(crv3.poles().begin(),crv3.poles().end(), crv4.poles().begin(),poles5.begin(),
    // [&](const auto &thick,const auto &cl){
    //     auto u_cl = gbs::extrema_curve_point(crv4,gbs::weight_projection(cl),1e-6).u;
    //     auto t_cl = crv4.value(u_cl,1); 
    //     std::array<double,3> n_cl = {-t_cl[1],t_cl[0],0};
    //     // auto n_cl = crv4.value(u_cl,2); 
    //     n_cl = n_cl / gbs::norm(n_cl);
    //     // auto ep = crv3.value(u_cl)[1];
    //     auto ep = thick[1];
    //     // auto thick_p = gbs::weight_projection(thick);
    //     // auto cl_p = gbs::weight_projection(cl);
    //     std::array<double,4> p{cl[0]+n_cl[0]*ep,cl[1]+n_cl[1]*ep,cl[2]+n_cl[2]*ep,thick[3]};
    //     return p;
    //     });

    gbs::BSCurveRational3d_d crv5(poles5,crv3.knotsFlats(),2);

    gbs::plot(crv3,crv4,crv5);

    occt_utils::to_iges(std::vector<Handle_Geom_Curve>{occt_utils::NURBSplineCurve(crv3),occt_utils::NURBSplineCurve(crv4)},"impeller_thickness.igs",1.);
    
}

TEST(tests_vtk_render, points)
{

    std::string line;
    std::ifstream myfile("../tests/in/e1098.dat");
    // std::ifstream myfile("../tests/in/e817.dat");
    if (myfile.is_open())
    {
        std::vector<std::array<double, 2>> pts;
        getline(myfile, line);
        while (getline(myfile, line))
        {
            std::istringstream iss(line);
            std::string::size_type sz; // alias of size_t

            double x = std::stod(line, &sz);
            double y = std::stod(line.substr(sz));
            pts.push_back({x, y});
        }
        myfile.close();

        auto crv = gbs::approx(pts, 5, gbs::KnotsCalcMode::CHORD_LENGTH,true);

        auto colors = vtkSmartPointer<vtkNamedColors>::New();
        auto pointActor = gbs::make_actor(pts,10.,true,colors->GetColor4d("Blue").GetData());
        gbs::plot(pointActor,crv);

    }
}

TEST(tests_vtk_render, surfNRUBS_points)
{
    std::vector<double> ku = {0.,0.,0.,1.,2.,3.,4.,4.,5.,5.,5.};
    std::vector<double> kv = {0.,0.,0.,1.,2.,3.,3.,3.};
    size_t p = 2;
    size_t q = 2;

    std::vector<std::array<double,4> > poles_t = {
        {0,2,4,1},{0,2,4,1}, {0,6,4,2},    {0,2,0,1},{0,2,0,1},
        {0,2,4,1},{0,2,4,1}, {0,6,4,2},    {0,2,0,1},{0,2,0,1},
        {0,2,4,1},{0,2,4,1}, {0,6,4,2},    {0,2,0,1},{0,2,0,1},
        {0,2,4,1},{4,6,8,2}, {12,24,12,6}, {4,6,0,2},{0,2,0,1},
        {0,2,4,1},{4,2,4,1}, {8,6,4,2},    {4,2,0,1},{0,2,0,1},
        {0,2,4,1},{4,2,4,1}, {8,6,4,2},    {4,2,0,1},{0,2,0,1},
        {0,2,4,1},{4,2,4,1}, {8,6,4,2},    {4,2,0,1},{0,2,0,1},
        {0,2,4,1},{4,2,4,1}, {8,6,4,2},    {4,2,0,1},{0,2,0,1}
                                                    };

    //Pij avec i inner loop
    std::vector<std::array<double,4> > poles(poles_t.size());
    int ni = 5 , nj =8;
    for (int i = 0; i < ni; i++)
    {
        for (int j = 0; j < nj; j++)
        {
            poles[j + nj * i] = poles_t[i + ni * j];
        }
    }

    gbs::BSSurfaceRational<double,3> srf(poles,ku,kv,p,q);

    auto pts = gbs::discretize(srf,20,30);

    auto colors = vtkSmartPointer<vtkNamedColors>::New();
    auto pointActor = gbs::make_actor(pts,5.,true,colors->GetColor4d("Blue").GetData());

    gbs::points_vector<double, 3> poles_;
    std::vector<double> weights;
    gbs::separate_weights(srf.poles(), poles_, weights);

    auto polesActor = gbs::make_actor(poles_,15.,true,colors->GetColor4d("Red").GetData());
    gbs::plot(pointActor,polesActor);


}

TEST(tests_vtk_render, surf_points)
{
    //Pij avec j inner loop
    // ---U--
    // |
    // V
    // |
    const gbs::points_vector<double,3> points =
    {
        {0,0,0},{1,0,0},
        {0,1,0},{1,1,1},
        {0,2,1},{2,1,0},
        {3,2,0},{3,2,0},
    };
    std::vector<double> ku = {0.,0.,1.,1.};
    std::vector<double> kv = {0.,0.,0.,0.5,1.,1.,1.};
    size_t p = 1;
    size_t q = 2;
    std::vector<double> u = {0.,1.};
    std::vector<double> v = {0.,0.33,0.66,1.};
    auto poles = gbs::build_poles(points,ku,kv,u,v,p,q);

    gbs::BSSurface<double,3> srf(poles,ku,kv,p,q) ;

    auto colors = vtkSmartPointer<vtkNamedColors>::New();

    gbs::plot(
        gbs::make_actor(points,10.,true,colors->GetColor4d("Blue").GetData()),
        gbs::make_actor(srf)
        );
}

#include <vtkActor.h>
#include <vtkAreaPicker.h>
#include <vtkCamera.h>
#include <vtkCellPicker.h>
#include <vtkDataSetMapper.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkExtractGeometry.h>
#include <vtkGlyph3D.h>
#include <vtkIdFilter.h>
#include <vtkIdTypeArray.h>
#include <vtkInteractorStyleTrackballActor.h>
#include <vtkNamedColors.h>
#include <vtkNew.h>
#include <vtkObjectFactory.h>
#include <vtkPlanes.h>
#include <vtkPointData.h>
#include <vtkPointSource.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkRendererCollection.h>
#include <vtkSmartPointer.h>
#include <vtkSphereSource.h>
#include <vtkPointSource.h>
#include <vtkUnstructuredGrid.h>
#include <vtkVertexGlyphFilter.h>
// /*
namespace
{
    template <typename T, size_t dim, bool rational>
    class InteractorStyleEditCurve : public vtkInteractorStyleTrackballActor
    {
    public:
        gbs::BSCurveGeneral<T,dim,rational> *p_crv_;
        vtkSmartPointer<vtkActor> ctr_polygon_lines; 
        vtkSmartPointer<vtkActor> crv_lines;
        
        vtkNew<vtkNamedColors> color;

        vtkPolyData* Data;
        vtkPolyData* GlyphData;

        vtkSmartPointer<vtkPolyDataMapper> MoveMapper;
        vtkSmartPointer<vtkActor> MoveActor;
        vtkSmartPointer<vtkSphereSource> MoveSphereSource;

        bool Move;
        vtkIdType SelectedPoint;
    public:
        static InteractorStyleEditCurve *New() // Can't make vtkStandardNewMacro works with templates
        {
            auto result = new InteractorStyleEditCurve<T,dim,rational>{};
            result->InitializeObjectBase();
            return result;
        }
        vtkTypeMacro(InteractorStyleEditCurve, vtkInteractorStyleTrackballActor);

        InteractorStyleEditCurve() : p_crv_{nullptr}
        {
            this->MoveSphereSource = vtkSmartPointer<vtkSphereSource>::New();
            this->MoveSphereSource->SetRadius(.1);
            this->MoveSphereSource->SetPhiResolution(36);
            this->MoveSphereSource->SetThetaResolution(36);
            this->MoveSphereSource->Update();

            this->MoveMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
            this->MoveMapper->SetInputConnection(
                this->MoveSphereSource->GetOutputPort());

            this->MoveActor = vtkSmartPointer<vtkActor>::New();
            this->MoveActor->SetMapper(this->MoveMapper);
            this->MoveActor->GetProperty()->SetOpacity(0.3);
            this->MoveActor->GetProperty()->SetColor(
                this->color->GetColor3d("Pink").GetData());
            // this->MoveActor->VisibilityOff();

            this->Move = false;
        }

        void OnMouseMove() override
        {
            if (!this->Move)
            {
                return;
            }
            vtkInteractorStyleTrackballActor::OnMouseMove();

            auto pos = this->MoveActor->GetPosition();
            gbs::point<T,dim> new_pole_pos{};
            for(auto i = 0; i < fmin(dim,3);i++)
            {
                new_pole_pos[i] = pos[i];
            }

            p_crv_->changePole(this->SelectedPoint,new_pole_pos);
            auto pts = gbs::discretize(*p_crv_, 30, 0.05);
            auto poles = p_crv_->poles();

            auto ctr_polygon_lines = gbs::make_polyline(poles,color->GetColor3d("Black").GetData());
            ctr_polygon_lines->GetProperty()->SetOpacity(0.3);
            ctr_polygon_lines->GetProperty()->SetLineWidth(3.f);
            auto crv_lines = gbs::make_polyline(pts,color->GetColor3d("Tomato").GetData());
            this->GetCurrentRenderer()->AddActor(ctr_polygon_lines);
            this->GetCurrentRenderer()->AddActor(crv_lines);
            this->GetCurrentRenderer()->RemoveActor(this->ctr_polygon_lines);
            this->GetCurrentRenderer()->RemoveActor(this->crv_lines);
            this->ctr_polygon_lines =ctr_polygon_lines;
            this->crv_lines = crv_lines;

            this->GetCurrentRenderer()->Render();
            this->GetCurrentRenderer()->GetRenderWindow()->Render();


        }

        void OnMiddleButtonUp() override
        {
            // Forward events
            vtkInteractorStyleTrackballActor::OnMiddleButtonUp();
            this->Move = false;
            this->MoveActor->VisibilityOff();

            this->Data->GetPoints()->SetPoint(this->SelectedPoint,
                                            this->MoveActor->GetPosition());
            this->Data->Modified();
            this->GetCurrentRenderer()->Render();
            this->GetCurrentRenderer()->GetRenderWindow()->Render();
        }

        void OnMiddleButtonDown() override
        {
            // Forward events
            vtkInteractorStyleTrackballActor::OnMiddleButtonDown();
            this->MoveActor->VisibilityOn();
            if (static_cast<vtkCellPicker*>(this->InteractionPicker)->GetPointId() >= 0)
            {
            vtkIdType id =
                dynamic_cast<vtkIdTypeArray*>(
                    this->GlyphData->GetPointData()->GetArray("InputPointIds"))
                    ->GetValue(static_cast<vtkCellPicker*>(this->InteractionPicker)
                                    ->GetPointId());
            // std::cout << "Id: " << id << std::endl;
            this->Move = true;
            this->SelectedPoint = id;

            double p[3];
            this->Data->GetPoint(id, p);
            // std::cout << "p: " << p[0] << " " << p[1] << " " << p[2] << std::endl;
            this->MoveActor->SetPosition(p);
            }

            this->GetCurrentRenderer()->AddActor(this->MoveActor);
            this->InteractionProp = this->MoveActor;
        }

        // void setCurve(gbs::BSCurveGeneral<T,dim,rational> *p_crv, vtkSmartPointer<vtkAssembly> p_actor)
        // {
        //     p_crv_ = p_crv;
        //     p_actor_ = p_actor;
        // }
    };
}
// */



TEST(tests_vtk_render, editCurve)
{
    std::vector<double> k = {0., 0., 0., 0., 1., 1., 1., 1.};
    std::vector<std::array<double,3> > poles =
    {
        {0.,1.,0.},
        {1.,2.,0.},
        {2.,2.,0.},
        {3.,0.,0.},
    };

    size_t p = 3;

    gbs::BSCurve3d_d crv(poles,k,p);


    // auto actor =gbs::make_actor(crv);
    auto pts = gbs::discretize(crv, 30, 0.01);
    // auto poles = crv.poles();
    vtkNew<vtkNamedColors> color;
    auto ctr_polygon_lines = gbs::make_polyline(poles,color->GetColor3d("Black").GetData());
    ctr_polygon_lines->GetProperty()->SetOpacity(0.3);
    auto crv_lines = gbs::make_polyline(pts,color->GetColor3d("Tomato").GetData());

    auto vtk_poles = gbs::make_vtkPoints(poles);
   

    vtkNew<vtkPolyData> input;
    input->SetPoints(vtk_poles);

    vtkNew<vtkSphereSource> glyphSource;
    // vtkNew<vtkPointSource> glyphSource;
    // glyphSource->SetRadius(2.5);

    glyphSource->SetRadius(0.1);
    glyphSource->SetThetaResolution(36);
    glyphSource->SetPhiResolution(36);
    glyphSource->Update();

    vtkNew<vtkGlyph3D> glyph3D;
    glyph3D->GeneratePointIdsOn();
    glyph3D->SetSourceConnection(glyphSource->GetOutputPort());
    glyph3D->SetInputData(input);
    glyph3D->SetScaleModeToDataScalingOff();
    // glyph3D->SetScaling(true);
    glyph3D->Update();

    // Create a mapper and actor
    vtkNew<vtkPolyDataMapper> mapper;
    mapper->SetInputConnection(glyph3D->GetOutputPort());

    vtkNew<vtkActor> actor;
    actor->SetMapper(mapper);
    actor->GetProperty()->SetColor(color->GetColor3d("Lime").GetData());
    actor->GetProperty()->SetOpacity(0.3);

    // Visualize
    vtkNew<vtkRenderer> renderer;
    vtkNew<vtkRenderWindow> renderWindow;
    renderWindow->AddRenderer(renderer);
    renderWindow->SetWindowName("MoveAGlyph");

    vtkNew<vtkRenderWindowInteractor> renderWindowInteractor;
    renderWindowInteractor->SetRenderWindow(renderWindow);

    renderer->AddActor(actor);
    renderer->AddActor(ctr_polygon_lines);
    renderer->AddActor(crv_lines);
    renderer->SetBackground(color->GetColor3d("White").GetData());

    renderWindow->Render();

    vtkNew<InteractorStyleEditCurve<double,3,false>> style;
    renderWindowInteractor->SetInteractorStyle(style);
    style->p_crv_ = &crv;
    style->Data = input;
    style->GlyphData = glyph3D->GetOutput();
    style->ctr_polygon_lines = ctr_polygon_lines;
    style->crv_lines = crv_lines;

    renderer->GetActiveCamera()->Zoom(0.9);
    renderWindowInteractor->Start();

}