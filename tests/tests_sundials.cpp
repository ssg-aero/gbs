#include <gtest/gtest.h>
#include "kinsol/kinsol.h" 
#include <nvector/nvector_serial.h> //N_Vector
#include <sundials/sundials_math.h> //NV_DATA_S
#include <sunlinsol/sunlinsol_spgmr.h>  //access to SPGMR SUNLinearSolver
#include <kinsol/kinsol_spils.h> // access to KINSpils interface

#include <gbslib/bscurve.h>
#include <gbslib/vecop.h>
#include <gbslib/bssinterp.h>

#include <algorithm>

using gbs::operator-;

TEST(tests_sundials, simple1)
{
    auto N = 2; // Taille du pb
    auto y0 = N_VNew_Serial(N);//Solution du Pb
    NV_DATA_S(y0)[0] = 2;
    NV_DATA_S(y0)[1] = 1;

    auto sc = N_VNew_Serial(N);//Le Pb ne requiert pas de mise à l'ech
    NV_DATA_S(sc)[0] = 1;
    NV_DATA_S(sc)[1] = 1;

    auto kin_mem = KINCreate();//allocation de la mem

//def de la solution à résoudre
    auto f = [](N_Vector u, N_Vector f_val, void *user_data) {
            auto p_u = N_VGetArrayPointer(u);     // pointer u vector data
            auto p_f = N_VGetArrayPointer(f_val); // pointer to udot vector data

            p_f[0] = p_u[0] * p_u[1] - 6.;
            p_f[1] = p_u[1] * sin(p_u[0])- 5.;

            return (0);
        };

    KINInit(
        kin_mem,
        f,
        y0);

    auto LS = SUNSPGMR(y0, 0, 0);
    KINSpilsSetLinearSolver(kin_mem, LS);

    /* Call KINSol and print output concentration profile */
    KINSol(kin_mem,        /* KINSol memory block */
                  y0,             /* initial guess on input; solution vector */
                  KIN_LINESEARCH, /* global strategy choice */
                  sc,             /* scaling vector for the variable cc */
                  sc);            /* scaling vector for function values fval */

    // Printing output.
    std::cout << "Final Value of y0 vector: \n";
    N_VPrint_Serial(y0);
    auto res = N_VNew_Serial(N);
    f(y0,res,nullptr);
    std::cout << "Final Residual: \n";
    N_VPrint_Serial(res);
}
TEST(tests_sundials, simple_proj)
{

    std::vector<double> k = {0., 0., 0., 1, 2, 3, 4, 5., 5., 5.};
    std::vector<std::array<double,3> > poles =
    {
        {0.,0.,0.},
        {0.,1.,0.},
        {1.,1.,0.},
        {1.,1.,1.},
        {1.,1.,2.},
        {3.,1.,1.},
        {0.,4.,1.},
    };
    size_t p = 2;
    double u = 2.3;

    auto c = gbs::BSCurve(poles,k,p);


    std::array<double,3> pt = c.value(u);
    auto delta=0.1;
    pt[1]+=delta;


    auto N = 1; // Taille du pb
    auto y0 = N_VNew_Serial(N);//Solution du Pb
    NV_DATA_S(y0)[0] = 2;

    auto sc = N_VNew_Serial(N);//Le Pb ne requiert pas de mise à l'ech
    NV_DATA_S(sc)[0] = 1;

    auto kin_mem = KINCreate();//allocation de la mem

    class UserData
    {
        public:
        UserData(const gbs::BSCurve<double,3> &crv,const std::array<double,3> &p) : c{crv}, pt{p} {}
        gbs::BSCurve<double, 3> c;
        std::array<double, 3> pt;
    };

    UserData data(c,pt);

    KINSetUserData(kin_mem, &data);

//def de la solution à résoudre
    auto f = [](N_Vector u, N_Vector f_val, void *user_data) {
            auto p_u = N_VGetArrayPointer(u);     // pointer u vector data
            auto p_f = N_VGetArrayPointer(f_val); // pointer to udot vector data
            // auto p_d = static_cast<UserData*>(user_data);
            auto p_d = (UserData*)(user_data);

            p_f[0] = gbs::norm( p_d->c.value(p_u[0]) - p_d->pt );

            return (0);
        };

    KINInit(
        kin_mem,
        f,
        y0);

    auto LS = SUNSPGMR(y0, 0, 0);
    KINSpilsSetLinearSolver(kin_mem, LS);

    
    KINSetErrFile(kin_mem, NULL);//As we look for minimum and only have root findig -> turn error msg off

    /* Call KINSol and print output concentration profile */
    KINSol(kin_mem,        /* KINSol memory block */
                  y0,             /* initial guess on input; solution vector */
                  KIN_LINESEARCH, /* global strategy choice */
                  sc,             /* scaling vector for the variable cc */
                  sc);            /* scaling vector for function values fval */

    // Printing output.
    std::cout << "Final Value of y0 vector: \n";
    N_VPrint_Serial(y0);
    auto res = N_VNew_Serial(N);
    f(y0,res,&data);
    std::cout << "Final Residual: \n";
    N_VPrint_Serial(res);

    ASSERT_NEAR(NV_DATA_S(res)[0],delta,1e-6);
}

template <typename T, size_t N>
class simple_solver
{
    void* kin_mem = nullptr;
    // void* user_data = nullptr;
    N_Vector y0;
    N_Vector sc;
    public:
    simple_solver()
    {
        kin_mem = KINCreate();
        y0 = N_VNew_Serial(N);
        sc = N_VNew_Serial(N);
        std::fill(
            N_VGetArrayPointer(sc),
            N_VGetArrayPointer(sc)+N,
            T(1.)
        );
    }
    template<typename Container>
    simple_solver(const Container &Y0) : simple_solver()
    {
        // kin_mem = KINCreate();
        // y0 = N_VNew_Serial(N);
        setInitVal(Y0);
    }
    ~simple_solver()
    {
        KINFree(&kin_mem);
        N_VDestroy(y0);
    }
    template<typename Container>
    auto setInitVal(const Container &Y0) -> void
    {
        if(Y0.size()!=N)
        {
            std::exception("wrong size");
        }
        std::copy(
            Y0.begin(),
            Y0.end(),
            N_VGetArrayPointer(y0)
            );
    }
    auto setUserData(void *data) -> void
    {
        KINSetUserData(kin_mem, data);
    }

    template <typename F>
    auto solve(F &f) -> std::array<T, N>
    {
        KINInit(
            kin_mem,
            f,
            y0);

        auto LS = SUNSPGMR(y0, 0, 0);
        KINSpilsSetLinearSolver(kin_mem, LS);

        // KINSetErrFile(kin_mem, NULL); //As we look for minimum and only have root findig -> turn error msg off

        /* Call KINSol and print output concentration profile */
        KINSol(kin_mem,        /* KINSol memory block */
               y0,             /* initial guess on input; solution vector */
               KIN_LINESEARCH, /* global strategy choice */
               sc,             /* scaling vector for the variable cc */
               sc);            /* scaling vector for function values fval */

        std::cout << "Final Value of y0 vector: \n";
        N_VPrint_Serial(y0);
        std::array<T, N> res;
        std::copy(
            N_VGetArrayPointer(y0),
            N_VGetArrayPointer(y0) + N,
            res.begin());
        SUNLinSolFree(LS);
        return res;
    }
};

TEST(tests_sundials, simple_projC_withsolver)
{

    std::vector<double> k = {0., 0., 0., 1, 2, 3, 4, 5., 5., 5.};
    std::vector<std::array<double,3> > poles =
    {
        {0.,0.,0.},
        {0.,1.,0.},
        {1.,1.,0.},
        {1.,1.,1.},
        {1.,1.,2.},
        {3.,1.,1.},
        {0.,4.,1.},
    };
    size_t p = 2;
    double u = 2.3;

    auto c = gbs::BSCurve(poles,k,p);


    std::array<double,3> pt = c.value(u);
    auto delta=0.1;
    pt[1]+=delta;

    class ExtreamPC
    {
    public:
        ExtreamPC(const gbs::BSCurve<double, 3> &crv, const std::array<double, 3> &p) : c{crv}, pt{p} {}
        gbs::BSCurve<double, 3> c;
        std::array<double, 3> pt;
    };

    ExtreamPC data(c, pt);

    simple_solver<double,1> solveur{};
    solveur.setInitVal(std::vector<double>({2.}));
    solveur.setUserData(&data);

    auto f = [](N_Vector u, N_Vector f_val, void *user_data) {
            auto p_u = N_VGetArrayPointer(u);     // pointer u vector data
            auto p_f = N_VGetArrayPointer(f_val); // pointer to udot vector data
            auto p_d = (ExtreamPC*)(user_data);

            p_f[0] = gbs::norm( p_d->c.value(p_u[0]) - p_d->pt );

            return (0);
        };
    
    auto sol = solveur.solve(f);
    ASSERT_NEAR(sol[0],u,1e-6);
    
}

TEST(tests_sundials, simple_projS_withsolver)
{

    const std::vector<std::array<double,3> > points =
    {
        {0,0,0},{1,0,0},
        {0,1,0},{1,1,1},
        {0,2,1},{2,1,0},
        {3,2,0},{3,2,0},
    };
    size_t p = 1;
    size_t q = 2;
    auto srf = gbs::interpolate(points,4,p,q,gbs::KnotsCalcMode::CHORD_LENGTH);
    auto u = 0.3;
    auto v = 0.7;


    std::array<double,3> pt = srf.value(u,v);
    auto delta=0.1;
    pt[1]+=delta;

    class ExtreamPS
    {
    public:
        ExtreamPS(const gbs::BSSurface<double, 3> &srf, const std::array<double, 3> &p) : s{srf}, pt{p} {}
        gbs::BSSurface<double, 3> s;
        std::array<double, 3> pt;
    };

    ExtreamPS data(srf, pt);

    simple_solver<double,2> solveur{};
    solveur.setInitVal(std::vector<double>({0.1,0.1}));
    solveur.setUserData(&data);

    auto f = [](N_Vector u, N_Vector f_val, void *user_data) {
            auto p_u = N_VGetArrayPointer(u);     // pointer u vector data
            auto p_f = N_VGetArrayPointer(f_val); // pointer to udot vector data
            auto p_d = (ExtreamPS*)(user_data);

            p_f[0] = gbs::norm( p_d->s.value(p_u[0],p_u[1]) - p_d->pt );
            p_f[1] = p_d->s.value(p_u[0],p_u[1])[0]-p_d->pt[0];
            // p_f[1] = gbs::norm( p_d->s.value(p_u[0],p_u[1]) - p_d->pt );

                    N_VPrint_Serial(u);
                    N_VPrint_Serial(f_val);

            return (0);
        };
    
    auto sol = solveur.solve(f);
    /****************************/
    /* Solver don't find minima */
    /* only roots               */
    /****************************/
    // ASSERT_NEAR(sol[0],u,1e-6);
    // ASSERT_NEAR(sol[1],v,1e-6);
    
}