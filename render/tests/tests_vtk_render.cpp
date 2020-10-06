#include <gtest/gtest.h>

#include <gbslib/bscanalysis.h>
#include <gbslib/bscbuild.h>
#include <gbslib/bscurve.h>
#include <gbslib/bscapprox.h>

#include <occt-utils/export.h>
#include <occt-utils/curvesbuild.h>

#include <render/vtkcurvesrender.h>

#include <vtkActor.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkNamedColors.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyLine.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSmartPointer.h>
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkTexture.h>
#include <vtkNew.h>
#include <vtkSphereSource.h>
#include <vtkPointSource.h>
#include <vtkVertexGlyphFilter.h>

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

    crv1.change_bounds(0,r1);
   
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
    crv3.change_bounds(0.,1.);

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

TEST(tests_vtk_render, dev)
{

    auto r1 = 3.;
    auto r2 = 1.;
    auto c = 25.;
    auto crv1 = gbs::build_ellipse<double, 3>(r1, r2, {r1, 0., 0.});
    crv1.trim(0.25, 0.5);
    crv1.reverse();

    crv1.change_bounds(0,r1);
   
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
    crv3.change_bounds(0.,1.);

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
     std::transform(std::next(crv3.poles().begin(),2),crv3.poles().end(),std::next(u.begin(),2),
     [&](const auto &pt_)
     {
         return gbs::extrema_PC(crv3,gbs::weight_projection( pt_ ),1e-6).u;
     });

    std::vector<std::array<double, 4>> poles5(crv3.poles().size());
    // std::transform(crv3.poles().begin(),crv3.poles().end(), crv4.poles().begin(),poles5.begin(),
    // [&](const auto &thick,const auto &cl){
    //     auto u_cl = gbs::extrema_PC(crv4,gbs::weight_projection(cl),1e-6).u;
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
        auto pointActor = gbs::make_actor(pts,5.,true,colors->GetColor3d("Lavender").GetData());
        gbs::plot(pointActor,crv);

    }
}