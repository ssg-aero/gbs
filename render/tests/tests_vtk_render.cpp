#include <gtest/gtest.h>

#include <gbslib/bscanalysis.h>
#include <gbslib/bscbuild.h>
#include <gbslib/bscurve.h>

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

    gbs::plot_curves(c_lsr);

}
    

}