#pragma once
#include <math_FunctionWithDerivative.hxx>
#include <math_FunctionSetWithDerivatives.hxx>
#include <math_MultipleVarFunctionWithHessian.hxx>
#include <math_Vector.hxx>
#include <math_Matrix.hxx>
namespace occt_utils
{
class Standard_EXPORT MultipleVarFunctionWithHessian : public math_MultipleVarFunctionWithHessian
{
public:
    MultipleVarFunctionWithHessian(double epsilon) {m_eps=epsilon;}
    bool 	Gradient(const math_Vector &X, math_Vector &G);
    bool 	Values  (const math_Vector &X, double &F, math_Vector &G);
    bool 	Values  (const math_Vector &X, double &F, math_Vector &G, math_Matrix &H);
    bool    Hessian (const math_Vector &X, math_Matrix &H);
    double eps() const;
    void setEps(double eps);
    double DX(const math_Vector &X,int i);
    double DX(const math_Vector &X,int i,int j);

protected:
    double m_eps;
    // bool debugMode;
    // int myGEvalcount;
    // int myHEvalcount;
    // int myEvalcount;
};
}