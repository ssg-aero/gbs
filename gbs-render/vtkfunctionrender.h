#include <vtkChartXY.h>
#include <vtkChartLegend.h>
#include <vtkContextScene.h>
#include <vtkContextView.h>
#include <vtkFloatArray.h>
#include <vtkNew.h>
#include <vtkPen.h>
#include <vtkPlot.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkTable.h>
#include <vector>
namespace gbs
{
    template< typename T>
    void plot_basis_funcs(size_t p, const std::vector<T> &u, const std::vector<std::vector<T>> &N_plot)
    {    
        size_t n =  N_plot.front().size();
        // Create a table with some points in it.
        vtkNew<vtkTable> table;
        vtkNew<vtkFloatArray> arrU;
        arrU->SetName("u");
        table->AddColumn(arrU);
        for(size_t i{}; i < n ; i++ )
        {
            vtkNew<vtkFloatArray> arr;
            arr->SetName( ("N"+std::to_string(i)+"_"+std::to_string(p)).c_str() );
            table->AddColumn(arr);    
        }
        // Fill in the table with some example values.
        int numPoints = N_plot.size();
        table->SetNumberOfRows(numPoints);
        for (int i = 0; i < numPoints; ++i)
        {
            table->SetValue(i, 0, u[i]);
            for(size_t j{}; j < n ; j++ )
            {
                table->SetValue(i, j + 1, N_plot[i][j]);
            }
        }

        // Set up the view
        vtkNew<vtkContextView> view;
        view->GetRenderWindow()->SetWindowName("Basis Functions");

        // Add multiple line plots, setting the colors etc.
        vtkNew<vtkChartXY> chart;
        chart->SetShowLegend(true);
        chart->GetLegend()->SetHorizontalAlignment(vtkChartLegend::LEFT);
        chart->GetLegend()->SetVerticalAlignment(vtkChartLegend::TOP);
        view->GetScene()->AddItem(chart);
        for(size_t j{1}; j <= n ; j++ )
        {
            auto line = chart->AddPlot(vtkChart::LINE);
            line->SetInputData(table, 0, j);
        }

        view->GetRenderWindow()->Render();
        view->GetInteractor()->Initialize();
        view->GetInteractor()->Start();
    }
}