#include "gbsBindCurves.h"

void gbs_bind_curves(py::module &m)
{
        gbs_bind_curves<double,1>(m);
        gbs_bind_curves<double,2>(m);
        gbs_bind_curves<double,3>(m);

        // Offsest without direction only have sense in 2d
        py::class_<gbs::CurveOffset<double,2,gbs::BSCfunction<double>>,std::shared_ptr<gbs::CurveOffset<double,2,gbs::BSCfunction<double>>>, gbs::Curve<double,2> >(m, "CurveOffset2d_bs")
        .def(py::init<const std::shared_ptr<gbs::Curve<double, 2>> &, const gbs::BSCfunction<double> &>())
        .def( "basisCurve", &gbs::CurveOffset< double, 2, gbs::BSCfunction<double> >::basisCurve )
        .def( "offset", &gbs::CurveOffset< double, 2, gbs::BSCfunction<double> >::offset )
        .def("__copy__", [](const gbs::CurveOffset<double,2,gbs::BSCfunction<double>> &self)
        { return gbs::CurveOffset<double,2,gbs::BSCfunction<double>>(self); })
        ;
        py::class_<gbs::CurveOffset<double,2,std::function<double(double,size_t)>>, std::shared_ptr<gbs::CurveOffset<double,2,std::function<double(double,size_t)>>>, gbs::Curve<double,2> >(m, "CurveOffset2d_func")
        .def(py::init<const std::shared_ptr<gbs::Curve<double, 2>> &, const std::function<double(double,size_t)> &>())
        .def( "basisCurve", &gbs::CurveOffset< double, 2, std::function<double(double,size_t)> >::basisCurve )
        .def( "offset", &gbs::CurveOffset< double, 2, std::function<double(double,size_t)> >::offset )
        .def("__copy__", [](const gbs::CurveOffset<double,2,std::function<double(double,size_t)>> &self)
        { return gbs::CurveOffset<double,2,std::function<double(double,size_t)>>(self); })
        ;
}