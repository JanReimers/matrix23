// File: matrix.cpp Unit tests for the Matrix<T> class.
#include "gtest/gtest.h"
#include <iostream>
#include <ranges>
#include "matrix23/matrix.hpp"

using std::cout;
using std::endl;
using matrix23::Vector;
using matrix23::Matrix;

class MatrixTests : public ::testing::Test
{
public:
    MatrixTests() = default;
    ~MatrixTests() override = default;

    template <std::ranges::range Range> static void print(Range v)
    {
        cout << "[";
        for (auto element : v) 
            std::cout << element << " ";

        cout << "]\n";
    }
    static auto constexpr print2D=[](auto rng){for(auto r:rng)print(r);};
    typedef std::initializer_list<double> il;
};


TEST_F(MatrixTests, Initialization)
{
    {
    Matrix<double, matrix23::FullRowMajorSubsciptor> A
    {{1, 2, 3},
     {4, 5, 6},
     {7, 8, 9}};
    EXPECT_EQ(A(0, 0), 1);
    EXPECT_EQ(A(0, 1), 2);
    EXPECT_EQ(A(0, 2), 3);
    EXPECT_EQ(A(1, 0), 4);
    EXPECT_EQ(A(1, 1), 5);
    EXPECT_EQ(A(1, 2), 6);
    EXPECT_EQ(A(2, 0), 7);
    EXPECT_EQ(A(2, 1), 8);
    EXPECT_EQ(A(2, 2), 9);
    EXPECT_EQ(A.row(0),il({1,2,3}));
    EXPECT_EQ(A.row(1),il({4,5,6}));
    EXPECT_EQ(A.row(2),il({7,8,9}));
    EXPECT_EQ(A.col(0),il({1,4,7}));
    EXPECT_EQ(A.col(1),il({2,5,8}));
    EXPECT_EQ(A.col(2),il({3,6,9}));
   
    }
    {
      Matrix<double, matrix23::FullColMajorSubsciptor> A
    {{1, 2, 3},
     {4, 5, 6},
     {7, 8, 9}};
    EXPECT_EQ(A(0, 0), 1);
    EXPECT_EQ(A(0, 1), 2);
    EXPECT_EQ(A(0, 2), 3);
    EXPECT_EQ(A(1, 0), 4);
    EXPECT_EQ(A(1, 1), 5);
    EXPECT_EQ(A(1, 2), 6);
    EXPECT_EQ(A(2, 0), 7);
    EXPECT_EQ(A(2, 1), 8);
    EXPECT_EQ(A(2, 2), 9);
    EXPECT_EQ(A.row(0),il({1,2,3}));
    EXPECT_EQ(A.row(1),il({4,5,6}));
    EXPECT_EQ(A.row(2),il({7,8,9}));
    EXPECT_EQ(A.col(0),il({1,4,7}));
    EXPECT_EQ(A.col(1),il({2,5,8}));
    EXPECT_EQ(A.col(2),il({3,6,9}));
    }
    {
    const Matrix<double, matrix23::UpperTriangularRowMajorSubsciptor> C
    {{1, 2, 3},
     {0, 5, 6},
     {0, 0, 9}};
    EXPECT_EQ(C(0, 0), 1);
    EXPECT_EQ(C(0, 1), 2);
    EXPECT_EQ(C(0, 2), 3);
    EXPECT_EQ(C(1, 0), 0); // Matrix has to be constant otherwise the non-const op(i,j) will fail.
    EXPECT_EQ(C(1, 1), 5);
    EXPECT_EQ(C(1, 2), 6);
    EXPECT_EQ(C(2, 0), 0); // Matrix has to be constant otherwise the non-const op(i,j) will fail.
    EXPECT_EQ(C(2, 1), 0); // Matrix has to be constant otherwise the non-const op(i,j) will fail.
    EXPECT_EQ(C(2, 2), 9);
    EXPECT_EQ(C.row(0),il({1,2,3}));
    EXPECT_EQ(C.row(1),il({5,6}));
    EXPECT_EQ(C.row(2),il({9}));
    EXPECT_EQ(C.col(0),il({1}));
    EXPECT_EQ(C.col(1),il({2,5}));
    EXPECT_EQ(C.col(2),il({3,6,9}));
    }
    {
        const Matrix<double, matrix23::UpperTriangularColMajorSubsciptor> C
        {{1, 2, 3},
        {0, 5, 6},
        {0, 0, 9}};
        EXPECT_EQ(C(0, 0), 1);
        EXPECT_EQ(C(0, 1), 2);
        EXPECT_EQ(C(0, 2), 3);
        EXPECT_EQ(C(1, 0), 0); // Matrix has to be constant otherwise the non-const op(i,j) will fail.
        EXPECT_EQ(C(1, 1), 5);
        EXPECT_EQ(C(1, 2), 6);
        EXPECT_EQ(C(2, 0), 0); // Matrix has to be constant otherwise the non-const op(i,j) will fail.
        EXPECT_EQ(C(2, 1), 0); // Matrix has to be constant otherwise the non-const op(i,j) will fail.
        EXPECT_EQ(C(2, 2), 9);
        EXPECT_EQ(C.row(0),il({1,2,3}));
        EXPECT_EQ(C.row(1),il({5,6}));
        EXPECT_EQ(C.row(2),il({9}));
        EXPECT_EQ(C.col(0),il({1}));
        EXPECT_EQ(C.col(1),il({2,5}));
        EXPECT_EQ(C.col(2),il({3,6,9}));
    }
    {
        const Matrix<double, matrix23::LowerTriangularColMajorSubsciptor> C
        {{1, 0, 0},
         {2, 4, 0},
         {3, 5, 6}};
        EXPECT_EQ(C(0, 0), 1);
        EXPECT_EQ(C(0, 1), 0);
        EXPECT_EQ(C(0, 2), 0);
        EXPECT_EQ(C(1, 0), 2); // Matrix has to be constant otherwise the non-const op(i,j) will fail.
        EXPECT_EQ(C(1, 1), 4);
        EXPECT_EQ(C(1, 2), 0);
        EXPECT_EQ(C(2, 0), 3); // Matrix has to be constant otherwise the non-const op(i,j) will fail.
        EXPECT_EQ(C(2, 1), 5); // Matrix has to be constant otherwise the non-const op(i,j) will fail.
        EXPECT_EQ(C(2, 2), 6);
        EXPECT_EQ(C.row(0),il({1}));
        EXPECT_EQ(C.row(1),il({2,4}));
        EXPECT_EQ(C.row(2),il({3,5,6}));
        EXPECT_EQ(C.col(0),il({1,2,3}));
        EXPECT_EQ(C.col(1),il({4,5}));
        EXPECT_EQ(C.col(2),il({6}));
    }
     {
        const Matrix<double, matrix23::LowerTriangularRowMajorSubsciptor> C
        {{1, 0, 0},
         {2, 4, 0},
         {3, 5, 6}};
        EXPECT_EQ(C(0, 0), 1);
        EXPECT_EQ(C(0, 1), 0);
        EXPECT_EQ(C(0, 2), 0);
        EXPECT_EQ(C(1, 0), 2); // Matrix has to be constant otherwise the non-const op(i,j) will fail.
        EXPECT_EQ(C(1, 1), 4);
        EXPECT_EQ(C(1, 2), 0);
        EXPECT_EQ(C(2, 0), 3); // Matrix has to be constant otherwise the non-const op(i,j) will fail.
        EXPECT_EQ(C(2, 1), 5); // Matrix has to be constant otherwise the non-const op(i,j) will fail.
        EXPECT_EQ(C(2, 2), 6);
        EXPECT_EQ(C.row(0),il({1}));
        EXPECT_EQ(C.row(1),il({2,4}));
        EXPECT_EQ(C.row(2),il({3,5,6}));
        EXPECT_EQ(C.col(0),il({1,2,3}));
        EXPECT_EQ(C.col(1),il({4,5}));
        EXPECT_EQ(C.col(2),il({6}));
    }
    
}

TEST_F(MatrixTests, FullMatrix)
{
    using matrix23::FullMatrix;
    FullMatrix mat(3, 4);
    mat(1, 2) = 5.0;
    EXPECT_EQ(mat(1, 2), 5.0);
    FullMatrix mat1{{1,2,3},
                    {4,5,6},
                    {7,8,9},
                    {10,11,12}};
    EXPECT_EQ(mat1(1, 2), 6.0);
    EXPECT_EQ(mat1(2, 1), 8.0);
    EXPECT_EQ(mat1(0, 0), 1.0);
    EXPECT_EQ(mat1(0, 1), 2.0);
    EXPECT_EQ(mat1(0, 2), 3.0);
    mat1.print();
    print(mat1.row(1)); // Should print: [4,5,6]
    print(mat1.col(2)); // Should print: [3,6,9,12]
    EXPECT_TRUE((mat1.row(1)==il{4,5,6}));
    EXPECT_TRUE((mat1.col(2)==il{3,6,9,12}));
}
TEST_F(MatrixTests, UpperTriangularMatrix)
{
    using matrix23::UpperTriangularMatrix;
    UpperTriangularMatrix mat(3, 4);
    mat(1, 2) = 5.0;
    EXPECT_EQ(mat(1, 2), 5.0);
    UpperTriangularMatrix mat1{{1,2,3},
                               {0,5,6},
                               {0,0,9},
                               {0,0,0}};
    EXPECT_EQ(mat1(1, 2), 6.0);
    EXPECT_EQ(mat1(0, 0), 1.0);
    EXPECT_EQ(mat1(0, 1), 2.0);
    EXPECT_EQ(mat1(0, 2), 3.0);
    mat1.print();
    print(mat1.row(1)); // Should print: [0,5,6]
    print(mat1.col(2)); // Should print: [3,6,9]
    EXPECT_TRUE((mat1.row(1)==il{5,6}));
    EXPECT_TRUE((mat1.col(2)==il{3,6,9}));
}
TEST_F(MatrixTests, DiagonalMatrix)
{
    using matrix23::DiagonalMatrix;
    DiagonalMatrix mat(3, 4);
    mat(1, 1) = 5.0;
    EXPECT_EQ(mat(1, 1), 5.0);
    DiagonalMatrix mat1{{1,0,0},
                        {0,5,0},
                        {0,0,9},
                        {0,0,0}};
    EXPECT_EQ(mat1(1, 1), 5.0);
    EXPECT_EQ(mat1(2, 2), 9.0);
    EXPECT_EQ(mat1(0, 0), 1.0);
    mat1.print();
    print(mat1.row(1)); // Should print: [0,5,0]
    print(mat1.col(2)); // Should print: [0,9]
    EXPECT_TRUE((mat1.row(1)==il{5}));
    EXPECT_TRUE((mat1.col(2)==il{9}));
}
TEST_F(MatrixTests, TriDiagonalMatrix)
{
    using matrix23::TriDiagonalMatrix;
    TriDiagonalMatrix mat(3, 3);
    mat(1, 0) = 5.0;
    mat(1, 1) = 6.0;
    mat(1, 2) = 7.0;
    EXPECT_EQ(mat(1, 0), 5.0);
    EXPECT_EQ(mat(1, 1), 6.0);
    EXPECT_EQ(mat(1, 2), 7.0);
    TriDiagonalMatrix mat1{{1,2,0},
                           {5,6,7},
                           {0,8,9}};
    EXPECT_EQ(mat1(1, 0), 5.0);
    EXPECT_EQ(mat1(1, 1), 6.0);
    EXPECT_EQ(mat1(1, 2), 7.0);
    EXPECT_EQ(mat1(2, 2), 9.0);
    EXPECT_EQ(mat1(0, 0), 1.0);
    mat1.print();
    print(mat1.row(1)); // Should print: [5,6,7]
    print(mat1.col(2)); // Should print: [7,9]
    EXPECT_TRUE((mat1.row(1)==il{5,6,7}));
    EXPECT_TRUE((mat1.col(2)==il{7,9}));
}

// #define RUN_SLOW_COMPILE_TESTS

#ifdef RUN_SLOW_COMPILE_TESTS
TEST_F(MatrixTests, DotProducts)
{
    using matrix23::Vector;
    using matrix23::TriDiagonalMatrix;
    Vector<double> v1{1,2,3,4,5};
    EXPECT_EQ(v1*v1,1*1+2*2+3*3+4*4+5*5);
    TriDiagonalMatrix A{{1,2,0,0,0},
                        {5,6,7,0,0},
                        {0,8,9,10,0,0},
                        {0,0,11,12,13},
                        {0,0,0,14,15}};
    EXPECT_EQ(A.row(0)*v1, 1*1 + 2*2); 
    EXPECT_EQ(A.row(1)*v1, 5*1 + 6*2 + 7*3); 
    EXPECT_EQ(A.row(2)*v1, 8*2 + 9*3 + 10*4); 
    EXPECT_EQ(A.row(3)*v1, 11*3 + 12*4 + 13*5); 
    EXPECT_EQ(A.row(4)*v1, 14*4 + 15*5); 
    EXPECT_EQ(A.col(0)*v1, 1*1 + 5*2); 
    EXPECT_EQ(A.col(1)*v1, 2*1 + 6*2 + 8*3); 
    EXPECT_EQ(A.col(2)*v1, 7*2 + 9*3 + 11*4); 
    EXPECT_EQ(A.col(3)*v1, 10*3 + 12*4 + 14*5); 
    EXPECT_EQ(A.col(4)*v1, 13*4 + 15*5); 
    
    EXPECT_EQ(v1*A.row(0), 1*1 + 2*2); 
    EXPECT_EQ(v1*A.row(1), 5*1 + 6*2 + 7*3); 
    EXPECT_EQ(v1*A.row(2), 8*2 + 9*3 + 10*4); 
    EXPECT_EQ(v1*A.row(3), 11*3 + 12*4 + 13*5); 
    EXPECT_EQ(v1*A.row(4), 14*4 + 15*5); 
    EXPECT_EQ(v1*A.col(0), 1*1 + 5*2); 
    EXPECT_EQ(v1*A.col(1), 2*1 + 6*2 + 8*3); 
    EXPECT_EQ(v1*A.col(2), 7*2 + 9*3 + 11*4); 
    EXPECT_EQ(v1*A.col(3), 10*3 + 12*4 + 14*5); 
    EXPECT_EQ(v1*A.col(4), 13*4 + 15*5);
 
    {
        Vector<double> Av(A * v1); // Matrix-vector multiplication
        EXPECT_EQ(Av(0), 1*1 + 2*2);
        EXPECT_EQ(Av(1), 5*1 + 6*2 + 7*3);
        EXPECT_EQ(Av(2), 8*2 + 9*3 + 10*4);
        EXPECT_EQ(Av(3), 11*3 + 12*4 + 13*5);
        EXPECT_EQ(Av(4), 14*4 + 15*5);
    }
    {
        Vector<double> Av=A * v1; // Matrix-vector multiplication
        EXPECT_EQ(Av(0), 1*1 + 2*2);
        EXPECT_EQ(Av(1), 5*1 + 6*2 + 7*3);
        EXPECT_EQ(Av(2), 8*2 + 9*3 + 10*4);
        EXPECT_EQ(Av(3), 11*3 + 12*4 + 13*5);
        EXPECT_EQ(Av(4), 14*4 + 15*5);
    }
    {
        Vector<double> vA = v1 * A; // Vector-Matrix multiplication
        EXPECT_EQ(vA(0), 1*1 + 5*2 + 0*3 + 0*4 + 0*5);
        EXPECT_EQ(vA(1), 2*1 + 6*2 + 8*3 + 0*4 + 0*5);
        EXPECT_EQ(vA(2), 0*1 + 7*2 + 9*3 + 11*4 + 0*5);
        EXPECT_EQ(vA(3), 0*1 + 0*2 + 10*3 + 12*4 + 14*5);
        EXPECT_EQ(vA(4), 0*1 + 0*2 + 0*3 + 13*4 + 15*5);
        print(vA); // Should print: [11,38,85,148,127]
    }
    EXPECT_EQ(v1*A*v1,11+2*38+3*85+4*148+5*127);
    
    
}

TEST_F(MatrixTests, MatriProducts)
{

    using::matrix23::TriDiagonalMatrix;
    using::matrix23::Vector;
    static_assert(matrix23::isMatrix<TriDiagonalMatrix>);

    TriDiagonalMatrix A{{1,2,0,0,0},
                        {5,6,7,0,0},
                        {0,8,9,10,0,0},
                        {0,0,11,12,13},
                        {0,0,0,14,15}};
    auto mpv=A*A; //MatriProductView
    EXPECT_EQ(mpv(0,0),11);
    EXPECT_EQ(mpv(0,1),14);
    EXPECT_EQ(mpv(0,2),14);
    EXPECT_EQ(mpv(1,0),35);
    EXPECT_EQ(mpv(1,1),102);
    EXPECT_EQ(mpv(1,2),105);
    EXPECT_EQ(mpv(1,3),70);
    EXPECT_EQ(mpv(2,0),40);
    EXPECT_EQ(mpv(2,1),120);
    EXPECT_EQ(mpv(2,2),247);
    EXPECT_EQ(mpv(2,3),210);
    EXPECT_EQ(mpv(2,4),130);
    EXPECT_EQ(mpv(3,1),88);
    EXPECT_EQ(mpv(3,2),231);
    EXPECT_EQ(mpv(3,3),436);
    EXPECT_EQ(mpv(3,4),351);
    EXPECT_EQ(mpv(4,2),154);
    EXPECT_EQ(mpv(4,3),378);
    EXPECT_EQ(mpv(4,4),407);

    auto mpvAAA=A*A*A;
    EXPECT_EQ(mpvAAA(0,0),81);

    Vector<double> v{1,2,3,4,5};
    Vector<double> vA=v*A;
    Vector<double> Av=A*v;
    print(vA);
    print(Av);
    typedef std::initializer_list<double> il;
    EXPECT_TRUE((vA==il{11, 38, 85, 148, 127}));
    EXPECT_TRUE((Av==il{5 , 38, 83, 146, 131}));
    EXPECT_EQ(v*A*v,1569);
    EXPECT_EQ(v*A*A*v,46799);
    EXPECT_EQ(v*A*A*A*v,1401625);
    EXPECT_EQ(v*A*A*A*A*v,42055783);
    EXPECT_EQ(v*A*A*A*A*A*v,1262119377); //This will slow down the compiler.
    
    using::matrix23::FullMatrix;
    FullMatrix An=A*A*A*A*A; //This will slow down the compiler.
}

TEST_F(MatrixTests, MatriOps)
{
    using::matrix23::TriDiagonalMatrix;
    using::matrix23::FullMatrix;
    using::matrix23::Vector;
    TriDiagonalMatrix A{{1,2,0,0,0},
                        {5,6,7,0,0},
                        {0,8,9,10,0,0},
                        {0,0,11,12,13},
                        {0,0,0,14,15}};

    TriDiagonalMatrix AA=A+A;
    AA.print();
    FullMatrix ApAA=A+A*A;
    ApAA.print();
    // print2D(plus.rows());
    FullMatrix A5=A*5,A6=6*A, A7=7.0*A, A8=A*8.0, A9=A/2, A10=A*A*9;
}
#endif