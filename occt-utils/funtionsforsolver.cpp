#include <occt-utils/funtionsforsolver.h>

namespace occt_utils
{
Standard_Boolean MultipleVarFunctionWithHessian::Gradient(const math_Vector &X, math_Vector &G){
    // myGEvalcount++;
    // if(debugMode) qDebug()<<"G Eval count :"<<myGEvalcount;
    // QVector<double> M;
    for(int i = 1 ; i <= X.Length() ; i++){
        G.Value(i) = DX(X,i);
        // if(debugMode) M<<G.Value(i);
        // std::cout << "G "<<i<<": "<<G.Value(i)<<std::endl;
    }
    // if(debugMode) qDebug()<<"G Values :"<<M;
    return true;
}

double MultipleVarFunctionWithHessian::DX(const math_Vector &X, int i) 
{

    math_Vector X1 = X;
    math_Vector X2 = X;
    X1.Value(i) = X.Value(i) - m_eps * 0.5;
    X2.Value(i) = X.Value(i) + m_eps * 0.5;
    double F1,F2;
    Value(X1,F1);
    Value(X2,F2);
    return (F2-F1) / m_eps;
}

double MultipleVarFunctionWithHessian::DX(const math_Vector &X, int i, int j){
    math_Vector X1 = X;
    math_Vector X2 = X;
    X1.Value(j) = X.Value(j) - m_eps * 0.5;
    X2.Value(j) = X.Value(j) + m_eps * 0.5;
    double F1,F2;
    F1 = DX(X1,i);
    F2 = DX(X2,i);
    return (F2-F1) / m_eps;
}

Standard_Boolean MultipleVarFunctionWithHessian::Hessian(const math_Vector &X, math_Matrix &H){
    // myHEvalcount++;
    // if(debugMode) qDebug()<<"H Eval count :"<<myHEvalcount;

    for(int i = 1 ; i <= X.Length() ; i++){
        // QVector<double> M;
        for(int j = 1 ; j <= X.Length() ; j++){
            H.Value(i,j) = DX(X,i,j);
            // if(debugMode) M<<H.Value(i,j);
            // std::cout << "H "<< i << "," << j  <<" : "<< H.Value(i,j) << std::endl; 
        }
        // if(debugMode) qDebug()<<"H["<<i<<"] Values :"<<M;
    }
    return true;
}

double MultipleVarFunctionWithHessian::eps() const
{
    return m_eps;
}

void MultipleVarFunctionWithHessian::setEps(double eps)
{
    m_eps = eps;
}

Standard_Boolean MultipleVarFunctionWithHessian::Values(const math_Vector &X, Standard_Real &F, math_Vector &G){
    return Value(X,F) && Gradient(X,G);
}

Standard_Boolean MultipleVarFunctionWithHessian::Values(const math_Vector &X, Standard_Real &F, math_Vector &G, math_Matrix &H){
    // myEvalcount++;
    // if(debugMode) qDebug()<<"Eval count :"<<myEvalcount;
    bool test=true;
    test = test &&  Values(X,F,G);
    // std::cout << 
    // QVector<double> M;
    // if(debugMode) {
        // for(int i = 1 ; i <= X.Length() ; i++){
    //         M<<X.Value(i);
            // std::cout << "X " << i << ": " << X.Value(i) << std::endl;
        // }
    //     if(debugMode) qDebug()<<"X Values :"<<M;
    //     qDebug()<<"F Value :"<<F;
    // }
    return test && Hessian(X,H);
}
}