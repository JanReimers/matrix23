// File: unittests/blas.cpp  Run some tests using the blas library.

#include "gtest/gtest.h"
#include <iostream>
#include "matrix23/matrix.hpp"
#include "matrix23/blas.hpp"

using std::cout;
using std::endl;
using matrix23::Vector;

class BlasTests : public ::testing::Test
{
public:
    BlasTests() = default;
    ~BlasTests() override = default;

    template <std::ranges::range Range> static void print(Range v)
    {
        cout << "[";
        for (auto element : v) 
            std::cout << element << " ";

        cout << "]\n";
    }
    static auto constexpr print2D=[](auto rng){for(auto r:rng)print(r);};
    typedef std::initializer_list<double> il;
    typedef std::initializer_list<il> ilil;
};

TEST_F(BlasTests,gemv)
{
    size_t nr=30,nc=40;
    {
        matrix23::FullMatrixCM<double> A(nr,nc,matrix23::random);
        Vector<double> x(nc,matrix23::random),y(nr);
        matrix23::gemv(1.0,A,x,0.0,y);
        EXPECT_EQ(A*x,y);
    }
    {
        matrix23::FullMatrixCM<double> A(nr,nc,matrix23::random);
        Vector<double> x(nc,matrix23::random);
        Vector<double> Ax=blasmv(A,x);
        EXPECT_EQ(A*x,Ax);
    }
    {
        matrix23::FullMatrixCM<double> A(nr,nc,matrix23::random);
        Vector<double> x(nr,matrix23::random);
        Vector<double> xA=blasvm(x,A);
        EXPECT_EQ(x*A,xA);
    }

    {
        matrix23::FullMatrixRM<double> A(nr,nc,matrix23::random);
        Vector<double> x(nc,matrix23::random),y(nr);
        matrix23::gemv(1.0,A,x,0.0,y);
        EXPECT_EQ(A*x,y);
    }
    {
        matrix23::FullMatrixRM<double> A(nr,nc,matrix23::random);
        Vector<double> x(nc,matrix23::random);
        Vector<double> Ax=blasmv(A,x);
        EXPECT_EQ(A*x,Ax);
    }
    {
        matrix23::FullMatrixRM<double> A(nr,nc,matrix23::random);
        Vector<double> x(nr,matrix23::random);
        Vector<double> xA=blasvm(x,A);
        EXPECT_EQ(x*A,xA);
    }

}

TEST_F(BlasTests,tpmv)
{
    size_t n=30;
    {
        matrix23::UpperTriangularMatrixCM<double> A(n,n,matrix23::random);
        Vector<double> x(n,matrix23::random),y(x);
        matrix23::tpmv(A,y);
        EXPECT_EQ(A*x,y);
    }
    {
        matrix23::UpperTriangularMatrixRM<double> A(n,n,matrix23::random);
        Vector<double> x(n,matrix23::random),y(x);
        matrix23::tpmv(A,y);
        EXPECT_EQ(A*x,y);
       
    }
    {
        matrix23::LowerTriangularMatrixCM<double> A(n,n,matrix23::random);
        Vector<double> x(n,matrix23::random),y(x);
        matrix23::tpmv(A,y);
        // EXPECT_EQ(A*x,y); flakey
        Vector<double> dx=A*x-y; 
        EXPECT_LT(sqrt(dx*dx),n*3e-16);
        // auto dx=A*x-y; //move problem somewhere.
    }
    {
        matrix23::LowerTriangularMatrixRM<double> A(n,n,matrix23::random);
        Vector<double> x(n,matrix23::random),y(x);
        matrix23::tpmv(A,y);
        EXPECT_EQ(A*x,y);
        Vector<double> dx=A*x-y; 
        EXPECT_LT(sqrt(dx*dx),n*2e-16);   
    }
}

TEST_F(BlasTests,gbmv)
{
    size_t n=30;
    for (size_t k=0;k<n-1;k++)
    {
        matrix23::SBandMatrix<double> A(n,k,matrix23::random);
        Vector<double> x(n,matrix23::random),y(n);
        matrix23::gbmv(1.0,A,x,0.0,y);
        EXPECT_EQ(A*x,y);
        Vector<double> Ax=matrix23::blasmv(A,x);
        EXPECT_EQ(A*x,Ax);

        Vector<double> xA=matrix23::blasvm(x,A);
        EXPECT_EQ(x*A,xA);
    }
}

TEST_F(BlasTests,ColMajor_gemm)
{
    using M=matrix23::FullMatrixCM<double>;
    size_t nr=30,k=35,nc=40;
    {
        M A(nr,k,matrix23::random);
        M B(k,nc,matrix23::random);
        M C(nr,nc);
        matrix23::gemm(1.0,A,B,0.0,C);
        // M delta=A*B-C;
        EXPECT_EQ(A*B,C);
    }
    {
        M A(nr,k,matrix23::random);
        M B(k,nc,matrix23::random);
        
        EXPECT_EQ(A*B,blasmm(A,B));
    }
}
// TEST_F(BlasTests,RowMajor_gemm)
// {
//     using M=matrix23::FullMatrixRM<double>;
//     // size_t nr=30,k=35,nc=40;
//     // size_t nr=1,k=4,nc=2;
//     {
//         // M A(nr,k)
//         // M B(k,nc)
//         // M A({{1,2,3,4}});
//         // M B({{4,1},{3,2},{2,3},{1,4}});
//         // M A({{1,2,3,4},
//         //      {4,3,2,1}});
//         // M B({{4},
//         //      {3},
//         //      {2},
//         //      {1}});
//         M A({{1,2},
//              {4,3}});
//         M B({{2},
//              {1}});
//         A.print();
//         B.print();
//         M C(A.nr(),B.nc());
//         matrix23::gemm(1.0,A,B,0.0,C);
//         C.print();

//         M C1=A*B;
//         C1.print();
//         // M delta=A*B-C;
//         // EXPECT_EQ(A*B,C);
//     }
//     // {
//     //     M A(nr,k,matrix23::one);
//     //     M B(k,nc,matrix23::one);
        
//     //     EXPECT_EQ(A*B,blasmm(A,B));
//     // }
// }

TEST_F(BlasTests,ColMajor_trmm)
{
    using F=matrix23::FullMatrixCM<double>;
    using UF=matrix23::UpperTriangularMatrixFCM<double>;
    using U=matrix23::UpperTriangularMatrixCM<double>;
    size_t nr=30,nc=40;
    {
        UF A=U(nr,nr,matrix23::random); //Do the fill on upper triagular packing and then copy into full packing.
        F B(nr,nc,matrix23::random);
        F AB=A*B; //Save the ranges answer before B gets modified.
        matrix23::trmm(1.0,A,B); //B=A*B
        EXPECT_EQ(AB,B);
    }
    {
        UF A=U(nc,nc,matrix23::random); //Do the fill on upper triagular packing and then copy into full packing.
        F B(nr,nc,matrix23::random);
        F BA=B*A; //Save the ranges answer before B gets modified.
        matrix23::trmm(1.0,B,A); //B=B*A
        F delta=BA-B;
        double d=fnorm(delta)/(nr*nc);
        EXPECT_LT(d,3e-16);
        cout << "d=" <<  d << endl;
    }
    using LF=matrix23::LowerTriangularMatrixFCM<double>;
    using L=matrix23::LowerTriangularMatrixCM<double>;
    {
        // nr=3;nc=4;
        LF A=L(nr,nr,matrix23::random); //Do the fill on lower triagular packing and then copy into full packing.
        F B(nr,nc,matrix23::random);
        F AB=A*B; //Save the ranges answer before B gets modified.
        matrix23::trmm(1.0,A,B); //B=A*B
        F delta=AB-B;
        double d=fnorm(delta)/(nr*nc);
        EXPECT_LT(d,3e-16);
        cout << "d=" <<  d << endl;
        // EXPECT_EQ(AB,B);
    }
    {
        LF A=L(nc,nc,matrix23::random); //Do the fill on lower triagular packing and then copy into full packing.
        F B(nr,nc,matrix23::random);
        F BA=B*A; //Save the ranges answer before B gets modified.
        matrix23::trmm(1.0,B,A); //B=A*B
        EXPECT_EQ(BA,B);
    }
   
}