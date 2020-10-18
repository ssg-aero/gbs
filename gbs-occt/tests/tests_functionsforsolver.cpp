#include <gtest/gtest.h>
#include <gbs-occt/funtionsforsolver.h>
#include <math_NewtonMinimum.hxx>

    class Standard_EXPORT EQ1 : public occt_utils::MultipleVarFunctionWithHessian
    {

    public:
        EQ1() : MultipleVarFunctionWithHessian(1.e-5) {}
        virtual Standard_Boolean Value(const math_Vector &X, Standard_Real &F) override
        {
            F = X(1)*X(1) + 2.*X(2)*X(2)+X(2);
            return true;
        }
        virtual Standard_Integer NbVariables() const override
        {
            return 2.;
        }
    };
TEST(tests_funtionsforsolver, multi_var_min)
{

        EQ1 eq;
        math_Vector X0(1,2);
        X0(1) = 1.;
        X0(2) = 1.;

        auto df12 = [](const math_Vector &X)
        {
            return 0.;
        };
        auto df21 = [](const math_Vector &X)
        {
            return 0.;
        };

        ASSERT_NEAR(eq.DX(X0,1),2.,1.e-5);
        ASSERT_NEAR(eq.DX(X0,2),5.,1.e-5);
        ASSERT_NEAR(eq.DX(X0,1,1),2.,1.e-5);
        ASSERT_NEAR(eq.DX(X0,2,2),4.,1.e-5);
        ASSERT_NEAR(eq.DX(X0,1,2),df12(X0),1.e-5);
        ASSERT_NEAR(eq.DX(X0,2,1),df12(X0),1.e-5);
        math_NewtonMinimum solver(eq,
                                Precision::Confusion(),
                                100,
                                1e-6,
                                true);


        solver.Perform(eq, X0);
        solver.Dump(std::cout);
}