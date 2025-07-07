// File: matrix__algebra.cpp Unit tests for the algebra with Matrix<T> classes.
#include "gtest/gtest.h"
#include <iostream>
#include <ranges>
#include "matrix23/matrix1.hpp"

using std::cout;
using std::endl;
using matrix23::Vector;

class MatrixAlgebraTests : public ::testing::Test
{
public:
    MatrixAlgebraTests() = default;
    ~MatrixAlgebraTests() override = default;

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


TEST_F(MatrixAlgebraTests, FullColMajor3x4)
{
    matrix23::FullMatrixCM<double> A({
        {1,2,3,4},
        {5,6,7,8},
        {9,10,11,12}});
    Vector<double> vr{1,2,3,4},vl{5,6,7};
    Vector<double> Av=A*vr;
    EXPECT_EQ(Av  ,(il{1*1+2*2+3*3+4*4,1*5+2*6+3*7+4*8,1*9+2*10+3*11+4*12}));
    EXPECT_EQ(A*vr,(il{1*1+2*2+3*3+4*4,1*5+2*6+3*7+4*8,1*9+2*10+3*11+4*12}));
    Vector<double> vA=vl*A;
    EXPECT_EQ(vA  ,(il{1*5+5*6+9*7,2*5+6*6+10*7,3*5+7*6+11*7,4*5+8*6+12*7}));
    EXPECT_EQ(vl*A,(il{1*5+5*6+9*7,2*5+6*6+10*7,3*5+7*6+11*7,4*5+8*6+12*7}));
    EXPECT_EQ(vl*A*vr,(1*5+5*6+9*7)*1+(2*5+6*6+10*7)*2+(3*5+7*6+11*7)*3+(4*5+8*6+12*7)*4);
}
TEST_F(MatrixAlgebraTests, FullRowMajor3x4)
{
    matrix23::FullMatrixRM<double> A({
        {1,2,3,4},
        {5,6,7,8},
        {9,10,11,12}});
    Vector<double> vr{1,2,3,4},vl{5,6,7};
    Vector<double> Av=A*vr;
    EXPECT_EQ(Av,(il{1*1+2*2+3*3+4*4,1*5+2*6+3*7+4*8,1*9+2*10+3*11+4*12}));
    EXPECT_EQ(A*vr,(il{1*1+2*2+3*3+4*4,1*5+2*6+3*7+4*8,1*9+2*10+3*11+4*12}));
    Vector<double> vA=vl*A;
    EXPECT_EQ(vA,(il{1*5+5*6+9*7,2*5+6*6+10*7,3*5+7*6+11*7,4*5+8*6+12*7}));
    EXPECT_EQ(vl*A,(il{1*5+5*6+9*7,2*5+6*6+10*7,3*5+7*6+11*7,4*5+8*6+12*7}));
    EXPECT_EQ(vl*A*vr,(1*5+5*6+9*7)*1+(2*5+6*6+10*7)*2+(3*5+7*6+11*7)*3+(4*5+8*6+12*7)*4);
}
TEST_F(MatrixAlgebraTests, FullColMajor4x3)
{
    matrix23::FullMatrixCM<double> A({
        {1,2,3},
        {4,5,6},
        {7,8,9},
        {10,11,12}});
    Vector<double> vl{1,2,3,4},vr{5,6,7};
    Vector<double> Av=A*vr;
    EXPECT_EQ(Av  ,(il{1*5+2*6+3*7,4*5+5*6+6*7,7*5+8*6+9*7,10*5+11*6+12*7}));
    EXPECT_EQ(A*vr,(il{1*5+2*6+3*7,4*5+5*6+6*7,7*5+8*6+9*7,10*5+11*6+12*7}));
    Vector<double> vA=vl*A;
    EXPECT_EQ(vA  ,(il{1*1+4*2+7*3+10*4,2*1+5*2+8*3+11*4,3*1+6*2+9*3+12*4}));
    EXPECT_EQ(vl*A,(il{1*1+4*2+7*3+10*4,2*1+5*2+8*3+11*4,3*1+6*2+9*3+12*4}));
    EXPECT_EQ(vl*A*vr,(1*1+4*2+7*3+10*4)*5 + (2*1+5*2+8*3+11*4)*6 + (3*1+6*2+9*3+12*4)*7);
}
TEST_F(MatrixAlgebraTests, FullRowMajor4x3)
{
    matrix23::FullMatrixRM<double> A({
        {1,2,3},
        {4,5,6},
        {7,8,9},
        {10,11,12}});
    Vector<double> vl{1,2,3,4},vr{5,6,7};
    Vector<double> Av=A*vr;
    EXPECT_EQ(Av,(il{1*5+2*6+3*7,4*5+5*6+6*7,7*5+8*6+9*7,10*5+11*6+12*7}));
    EXPECT_EQ(A*vr,(il{1*5+2*6+3*7,4*5+5*6+6*7,7*5+8*6+9*7,10*5+11*6+12*7}));
    Vector<double> vA=vl*A;
    EXPECT_EQ(vA,(il{1*1+4*2+7*3+10*4,2*1+5*2+8*3+11*4,3*1+6*2+9*3+12*4}));
    EXPECT_EQ(vl*A,(il{1*1+4*2+7*3+10*4,2*1+5*2+8*3+11*4,3*1+6*2+9*3+12*4}));
    EXPECT_EQ(vl*A*vr,(1*1+4*2+7*3+10*4)*5 + (2*1+5*2+8*3+11*4)*6 + (3*1+6*2+9*3+12*4)*7);
}
TEST_F(MatrixAlgebraTests, UpperTriangularColMajor3x4)
{
    matrix23::UpperTriangularMatrixCM<double> A({
        {1,2,3,4},
        {0,6,7,8},
        {0,0,11,12}});
    Vector<double> vr{1,2,3,4},vl{5,6,7};
    Vector<double> Av=A*vr;
    EXPECT_EQ(Av  ,(il{1*1+2*2+3*3+4*4,1*0+2*6+3*7+4*8,1*0+2*0+3*11+4*12}));
    EXPECT_EQ(A*vr,(il{1*1+2*2+3*3+4*4,1*0+2*6+3*7+4*8,1*0+2*0+3*11+4*12}));
    Vector<double> vA=vl*A;
    EXPECT_EQ(vA  ,(il{1*5+0*6+0*7,2*5+6*6+0*7,3*5+7*6+11*7,4*5+8*6+12*7}));
    EXPECT_EQ(vl*A,(il{1*5+0*6+0*7,2*5+6*6+0*7,3*5+7*6+11*7,4*5+8*6+12*7}));
    EXPECT_EQ(vl*A*vr,(1*5+0*6+0*7)*1+(2*5+6*6+0*7)*2+(3*5+7*6+11*7)*3+(4*5+8*6+12*7)*4);
}
TEST_F(MatrixAlgebraTests, UpperTriangularRowMajor3x4)
{
    matrix23::UpperTriangularMatrixRM<double> A({
        {1,2,3,4},
        {0,6,7,8},
        {0,0,11,12}});
    Vector<double> vr{1,2,3,4},vl{5,6,7};
    Vector<double> Av=A*vr;
    EXPECT_EQ(Av  ,(il{1*1+2*2+3*3+4*4,1*0+2*6+3*7+4*8,1*0+2*0+3*11+4*12}));
    EXPECT_EQ(A*vr,(il{1*1+2*2+3*3+4*4,1*0+2*6+3*7+4*8,1*0+2*0+3*11+4*12}));
    Vector<double> vA=vl*A;
    EXPECT_EQ(vA  ,(il{1*5+0*6+0*7,2*5+6*6+0*7,3*5+7*6+11*7,4*5+8*6+12*7}));
    EXPECT_EQ(vl*A,(il{1*5+0*6+0*7,2*5+6*6+0*7,3*5+7*6+11*7,4*5+8*6+12*7}));
    EXPECT_EQ(vl*A*vr,(1*5+0*6+0*7)*1+(2*5+6*6+0*7)*2+(3*5+7*6+11*7)*3+(4*5+8*6+12*7)*4);
}
TEST_F(MatrixAlgebraTests, UpperTriangularColMajor4x3)
{
    matrix23::UpperTriangularMatrixCM<double> A({
        {1,2,3},
        {0,5,6},
        {0,0,9},
        {0,0,0}});
    Vector<double> vl{1,2,3,4},vr{5,6,7};
    Vector<double> Av=A*vr;
    EXPECT_EQ(Av  ,(il{1*5+2*6+3*7,0*5+5*6+6*7,0*5+0*6+9*7,0*5+0*6+0*7}));
    EXPECT_EQ(A*vr,(il{1*5+2*6+3*7,0*5+5*6+6*7,0*5+0*6+9*7,0*5+0*6+0*7}));
    Vector<double> vA=vl*A;
    EXPECT_EQ(vA  ,(il{1*1+0*2+0*3+0*4,2*1+5*2+0*3+0*4,3*1+6*2+9*3+0*4}));
    EXPECT_EQ(vl*A,(il{1*1+0*2+0*3+0*4,2*1+5*2+0*3+0*4,3*1+6*2+9*3+0*4}));
    EXPECT_EQ(vl*A*vr,(1*1+0*2+0*3+0*4)*5 + (2*1+5*2+0*3+0*4)*6 + (3*1+6*2+9*3+0*4)*7);
}
TEST_F(MatrixAlgebraTests, UpperTriangularRowMajor4x3)
{
    matrix23::UpperTriangularMatrixRM<double> A({
        {1,2,3},
        {0,5,6},
        {0,0,9},
        {0,0,0}});
    Vector<double> vl{1,2,3,4},vr{5,6,7};
    Vector<double> Av=A*vr;
    EXPECT_EQ(Av  ,(il{1*5+2*6+3*7,0*5+5*6+6*7,0*5+0*6+9*7,0*5+0*6+0*7}));
    EXPECT_EQ(A*vr,(il{1*5+2*6+3*7,0*5+5*6+6*7,0*5+0*6+9*7,0*5+0*6+0*7}));
    Vector<double> vA=vl*A;
    EXPECT_EQ(vA  ,(il{1*1+0*2+0*3+0*4,2*1+5*2+0*3+0*4,3*1+6*2+9*3+0*4}));
    EXPECT_EQ(vl*A,(il{1*1+0*2+0*3+0*4,2*1+5*2+0*3+0*4,3*1+6*2+9*3+0*4}));
    EXPECT_EQ(vl*A*vr,(1*1+0*2+0*3+0*4)*5 + (2*1+5*2+0*3+0*4)*6 + (3*1+6*2+9*3+0*4)*7);
}
TEST_F(MatrixAlgebraTests, LowerTriangularColMajor3x4)
{
    matrix23::LowerTriangularMatrixCM<double> A({
        {1,0 ,0 ,0},
        {5,6 ,0 ,0},
        {9,10,11,0}});
    Vector<double> vr{1,2,3,4},vl{5,6,7};
    Vector<double> Av=A*vr;
    EXPECT_EQ(Av  ,(il{1*1+2*0+3*0+4*0,1*5+2*6+3*0+4*0,1*9+2*10+3*11+4*0}));
    EXPECT_EQ(A*vr,(il{1*1+2*0+3*0+4*0,1*5+2*6+3*0+4*0,1*9+2*10+3*11+4*0}));
    Vector<double> vA=vl*A;
    EXPECT_EQ(vA  ,(il{1*5+5*6+9*7,0*5+6*6+10*7,0*5+0*6+11*7,0*5+0*6+0*7}));
    EXPECT_EQ(vl*A,(il{1*5+5*6+9*7,0*5+6*6+10*7,0*5+0*6+11*7,0*5+0*6+0*7}));
    EXPECT_EQ(vl*A*vr,(1*5+5*6+9*7)*1+(0*5+6*6+10*7)*2+(0*5+0*6+11*7)*3+(0*5+0*6+0*7)*4);
}
TEST_F(MatrixAlgebraTests, LowerTriangularRowMajor3x4)
{
    matrix23::LowerTriangularMatrixRM<double> A({
        {1,0 ,0 ,0},
        {5,6 ,0 ,0},
        {9,10,11,0}});
    Vector<double> vr{1,2,3,4},vl{5,6,7};
    Vector<double> Av=A*vr;
    EXPECT_EQ(Av  ,(il{1*1+2*0+3*0+4*0,1*5+2*6+3*0+4*0,1*9+2*10+3*11+4*0}));
    EXPECT_EQ(A*vr,(il{1*1+2*0+3*0+4*0,1*5+2*6+3*0+4*0,1*9+2*10+3*11+4*0}));
    Vector<double> vA=vl*A;
    EXPECT_EQ(vA  ,(il{1*5+5*6+9*7,0*5+6*6+10*7,0*5+0*6+11*7,0*5+0*6+0*7}));
    EXPECT_EQ(vl*A,(il{1*5+5*6+9*7,0*5+6*6+10*7,0*5+0*6+11*7,0*5+0*6+0*7}));
    EXPECT_EQ(vl*A*vr,(1*5+5*6+9*7)*1+(0*5+6*6+10*7)*2+(0*5+0*6+11*7)*3+(0*5+0*6+0*7)*4);
}
TEST_F(MatrixAlgebraTests, LowerTriangularColMajor4x3)
{
    matrix23::LowerTriangularMatrixCM<double> A({
        {1,0,0},
        {4,5,0},
        {7,8,9},
        {10,11,12}});
    Vector<double> vl{1,2,3,4},vr{5,6,7};
    Vector<double> Av=A*vr;
    EXPECT_EQ(Av  ,(il{1*5+0*6+0*7,4*5+5*6+6*0,7*5+8*6+9*7,10*5+11*6+12*7}));
    EXPECT_EQ(A*vr,(il{1*5+0*6+0*7,4*5+5*6+6*0,7*5+8*6+9*7,10*5+11*6+12*7}));
    Vector<double> vA=vl*A;
    EXPECT_EQ(vA  ,(il{1*1+4*2+7*3+10*4,0*1+5*2+8*3+11*4,0*1+0*2+9*3+12*4}));
    EXPECT_EQ(vl*A,(il{1*1+4*2+7*3+10*4,0*1+5*2+8*3+11*4,0*1+0*2+9*3+12*4}));
    EXPECT_EQ(vl*A*vr,(1*1+4*2+7*3+10*4)*5 + (0*1+5*2+8*3+11*4)*6 + (0*1+0*2+9*3+12*4)*7);
}
TEST_F(MatrixAlgebraTests, LowerTriangularRowMajor4x3)
{
    matrix23::LowerTriangularMatrixRM<double> A({
        {1,0,0},
        {4,5,0},
        {7,8,9},
        {10,11,12}});
    Vector<double> vl{1,2,3,4},vr{5,6,7};
    Vector<double> Av=A*vr;
    EXPECT_EQ(Av  ,(il{1*5+0*6+0*7,4*5+5*6+6*0,7*5+8*6+9*7,10*5+11*6+12*7}));
    EXPECT_EQ(A*vr,(il{1*5+0*6+0*7,4*5+5*6+6*0,7*5+8*6+9*7,10*5+11*6+12*7}));
    Vector<double> vA=vl*A;
    EXPECT_EQ(vA  ,(il{1*1+4*2+7*3+10*4,0*1+5*2+8*3+11*4,0*1+0*2+9*3+12*4}));
    EXPECT_EQ(vl*A,(il{1*1+4*2+7*3+10*4,0*1+5*2+8*3+11*4,0*1+0*2+9*3+12*4}));
    EXPECT_EQ(vl*A*vr,(1*1+4*2+7*3+10*4)*5 + (0*1+5*2+8*3+11*4)*6 + (0*1+0*2+9*3+12*4)*7);
}
TEST_F(MatrixAlgebraTests, Diagonal3x4)
{
    matrix23::DiagonalMatrix<double> A({
        {1,0,0 ,0},
        {0,6,0 ,0},
        {0,0,11,0}});
    Vector<double> vr{1,2,3,4},vl{5,6,7};
    Vector<double> Av=A*vr;
    EXPECT_EQ(Av  ,(il{1*1,2*6,3*11}));
    EXPECT_EQ(A*vr,(il{1*1,2*6,3*11}));
    Vector<double> vA=vl*A;
    EXPECT_EQ(vA  ,(il{1*5,6*6,7*11,0}));
    EXPECT_EQ(vl*A,(il{1*5,6*6,7*11,0}));
    EXPECT_EQ(vl*A*vr,1*5*1 + 6*6*2 + 7*11*3 + 0);
}
TEST_F(MatrixAlgebraTests, Diagonal4x3)
{
    matrix23::DiagonalMatrix<double> A({
        {1,0,0},
        {0,5,0},
        {0,0,9},
        {0,0,0}});
    Vector<double> vl{1,2,3,4},vr{5,6,7};
    Vector<double> Av=A*vr;
    EXPECT_EQ(Av  ,(il{1*5,5*6,9*7,0}));
    EXPECT_EQ(A*vr,(il{1*5,5*6,9*7,0}));
    Vector<double> vA=vl*A;
    EXPECT_EQ(vA  ,(il{1*1,2*5,3*9}));
    EXPECT_EQ(vl*A,(il{1*1,2*5,3*9}));
    EXPECT_EQ(vl*A*vr,1*1*5 + 2*5*6 +3*9*7 + 0);
}
TEST_F(MatrixAlgebraTests, SBandk0)
{
    size_t k=0;
    matrix23::SBandMatrix<double> A({
        {1,0,0,0,0,0},
        {0,2,0,0,0,0},
        {0,0,3,0,0,0},
        {0,0,0,4,0,0},
        {0,0,0,0,5,0},
        {0,0,0,0,0,6}},k);
    Vector<double> v{7,8,9,10,11,12};
    EXPECT_EQ(A*v,(il{1*7,2*8,3*9,4*10,5*11,6*12}));
    EXPECT_EQ(v*A,(il{1*7,2*8,3*9,4*10,5*11,6*12}));
    EXPECT_EQ(v*A*v,1*7*7+2*8*8+3*9*9+4*10*10+5*11*11+6*12*12);
}
TEST_F(MatrixAlgebraTests, SBandk1)
{
    size_t k=1;
    matrix23::SBandMatrix<double> A({
        {1 , 7, 0, 0, 0, 0},
        {12, 2, 8, 0, 0, 0},
        {0 ,13, 3, 9, 0, 0},
        {0 , 0,14, 4,10, 0},
        {0 , 0, 0,15, 5,11},
        {0 , 0, 0, 0,16, 6}},k);
    Vector<double> v{7,8,9,10,11,12};
    EXPECT_EQ(A*v,(il{1*7+7*8,12*7+2*8+8*9,13*8+3*9+9*10,14*9+4*10+10*11,15*10+5*11+11*12,16*11+6*12}));
    EXPECT_EQ(v*A,(il{1*7+12*8,7*7+2*8+13*9,8*8+3*9+14*10,9*9+4*10+15*11,10*10+5*11+16*12,11*11+6*12}));
    EXPECT_EQ(v*A*v,(1*7+7*8)*7+(12*7+2*8+8*9)*8+(13*8+3*9+9*10)*9+(14*9+4*10+10*11)*10+(15*10+5*11+11*12)*11+(16*11+6*12)*12);
}
TEST_F(MatrixAlgebraTests, SBandk2)
{
    size_t k=2;
    matrix23::SBandMatrix<double> A({
        {1 , 7,17, 0, 0, 0},
        {12, 2, 8,18, 0, 0},
        {21,13, 3, 9,19, 0},
        {0 ,22,14, 4,10,20},
        {0 , 0,23,15, 5,11},
        {0 , 0, 0,24,16, 6}},k);
    Vector<double> v{7,8,9,10,11,12};
    EXPECT_EQ(A*v,(il{1*7+7*8+17*9,12*7+2*8+8*9+18*10,21*7+13*8+3*9+9*10+19*11,22*8+14*9+4*10+10*11+20*12,23*9+15*10+5*11+11*12,24*10+16*11+6*12}));
    EXPECT_EQ(v*A,(il{1*7+12*8+21*9,7*7+2*8+13*9+22*10,17*7+8*8+3*9+14*10+23*11,18*8+9*9+4*10+15*11+24*12,19*9+10*10+5*11+16*12,20*10+11*11+6*12}));
    EXPECT_EQ(v*A*v,(1*7+7*8+17*9)*7+(12*7+2*8+8*9+18*10)*8+(21*7+13*8+3*9+9*10+19*11)*9+(22*8+14*9+4*10+10*11+20*12)*10+(23*9+15*10+5*11+11*12)*11+(24*10+16*11+6*12)*12);
}
TEST_F(MatrixAlgebraTests, SBandk3)
{
    size_t k=3;
    matrix23::SBandMatrix<double> A({
        {1 , 7,17,25, 0, 0},
        {12, 2, 8,18,26, 0},
        {21,13, 3, 9,19,27},
        {28,22,14, 4,10,20},
        {0 ,29,23,15, 5,11},
        {0 , 0,30,24,16, 6}},k);
    Vector<double> v{7,8,9,10,11,12};
    EXPECT_EQ(A*v,(il{466, 638, 901, 888, 776, 758}));
    EXPECT_EQ(v*A,(il{572, 721, 963, 893, 726, 636}));
    EXPECT_EQ(v*A*v,42987);
}
TEST_F(MatrixAlgebraTests, SBandk4)
{
    size_t k=4;
    matrix23::SBandMatrix<double> A({
        {1 , 7,17,25,31, 0},
        {12, 2, 8,18,26,32},
        {21,13, 3, 9,19,27},
        {28,22,14, 4,10,20},
        {33,29,23,15, 5,11},
        {0 ,34,30,24,16, 6}},k);
    Vector<double> v{7,8,9,10,11,12};
    EXPECT_EQ(A*v,(il{807,1022,901,888,1007,1030}));
    EXPECT_EQ(v*A,(il{935,1129,963,893,943,892}));
    EXPECT_EQ(v*A*v,54251);
}
TEST_F(MatrixAlgebraTests, SBandk5)
{
    size_t k=5;
    matrix23::SBandMatrix<double> A({
        {1 , 7,17,25,31,35},
        {12, 2, 8,18,26,32},
        {21,13, 3, 9,19,27},
        {28,22,14, 4,10,20},
        {33,29,23,15, 5,11},
        {36,34,30,24,16, 6}},k);
    Vector<double> v{7,8,9,10,11,12};
    EXPECT_EQ(A*v,(il{1227,1022,901,888,1007,1282 }));
    EXPECT_EQ(v*A,(il{1367,1129,963,893,943,1137}));
    EXPECT_EQ(v*A*v,60215);
}
TEST_F(MatrixAlgebraTests, MatrixMultiplyFF)
{
    {
        matrix23::FullMatrixCM<double> A({
            {1,2,3,4},
            {5,6,7,8},
            {9,10,11,12}});
        matrix23::FullMatrixCM<double> B({
            {1,2,3},
            {4,5,6},
            {7,8,9},
            {10,11,12}});
        matrix23::FullMatrixCM<double> C=A*B;
        EXPECT_EQ(C.row(0),(il{70,80,90}));
        EXPECT_EQ(C.row(1),(il{158,184,210}));
        EXPECT_EQ(C.row(2),(il{246,288,330}));
        //matrix23::DiagonalMatrix<double> D=A*B; //Compile error
    }
    {
        matrix23::FullMatrixRM<double> A({
            {1,2,3,4},
            {5,6,7,8},
            {9,10,11,12}});
        matrix23::FullMatrixRM<double> B({
            {1,2,3},
            {4,5,6},
            {7,8,9},
            {10,11,12}});
        matrix23::FullMatrixRM<double> C=A*B;
        EXPECT_EQ(C.row(0),(il{70,80,90}));
        EXPECT_EQ(C.row(1),(il{158,184,210}));
        EXPECT_EQ(C.row(2),(il{246,288,330}));
    }
}
TEST_F(MatrixAlgebraTests, MatrixMultiplyFL)
{
    matrix23::FullMatrixCM<double> A({
        {1,2,3,4},
        {5,6,7,8},
        {9,10,11,12}});
    matrix23::LowerTriangularMatrixCM<double> B({
        {1,0,0},
        {4,5,0},
        {7,8,9},
        {10,11,12}});
    matrix23::FullMatrixCM<double> C=A*B;
    EXPECT_EQ(C.row(0),(il{70,78,75}));
    EXPECT_EQ(C.row(1),(il{158,174,159}));
    EXPECT_EQ(C.row(2),(il{246,270,243}));
}
TEST_F(MatrixAlgebraTests, MatrixMultiplyUF)
{
   matrix23::UpperTriangularMatrixCM<double> A({
        {1,2,3,4},
        {0,6,7,8},
        {0,0,11,12}});
    matrix23::FullMatrixCM<double> B({
        {1,2,3},
        {4,5,6},
        {7,8,9},
        {10,11,12}});
    matrix23::FullMatrixCM<double> C=A*B;
    EXPECT_EQ(C.row(0),(il{70,80,90}));
    EXPECT_EQ(C.row(1),(il{153,174,195}));
    EXPECT_EQ(C.row(2),(il{197,220,243}));
}
TEST_F(MatrixAlgebraTests, MatrixMultiplyDF)
{
    matrix23::DiagonalMatrix<double> A({
        {1,0,0,0},
        {0,6,0,0},
        {0,0,11,0}});
    matrix23::FullMatrixCM<double> B({
        {1,2,3},
        {4,5,6},
        {7,8,9},
        {10,11,12}});
    matrix23::FullMatrixCM<double> C=A*B;
    EXPECT_EQ(C.row(0),(il{1,2,3}));
    EXPECT_EQ(C.row(1),(il{24,30,36}));
    EXPECT_EQ(C.row(2),(il{77,88,99}));
}
TEST_F(MatrixAlgebraTests, MatrixMultiplyLL)
{
    matrix23::LowerTriangularMatrixCM<double> A({
        {1,0,0,0},
        {5,6,0,0},
        {9,10,11,0}});
    matrix23::LowerTriangularMatrixCM<double> B({
        {1,0,0},
        {4,5,0},
        {7,8,9},
        {10,11,12}});
    matrix23::LowerTriangularMatrixCM<double> C=A*B;
    EXPECT_EQ(C.row(0),(il{1}));
    EXPECT_EQ(C.row(1),(il{29,30}));
    EXPECT_EQ(C.row(2),(il{126,138,99}));
}
TEST_F(MatrixAlgebraTests, MatrixMultiplyUU)
{
    matrix23::UpperTriangularMatrixCM<double> A({
        {1,2,3,4},
        {0,6,7,8},
        {0,0,11,12}});
    matrix23::UpperTriangularMatrixCM<double> B({
        {1,2,3},
        {0,5,6},
        {0,0,9},
        {0,0,0}});
    matrix23::UpperTriangularMatrixCM<double> C=A*B;
    EXPECT_EQ(C.row(0),(il{1,12,42}));
    EXPECT_EQ(C.row(1),(il{30,99}));
    EXPECT_EQ(C.row(2),(il{99}));
  
}
TEST_F(MatrixAlgebraTests, MatrixMultiplyDL)
{
    matrix23::DiagonalMatrix<double> A({
        {1,0,0,0},
        {0,6,0,0},
        {0,0,11,0}});
    matrix23::LowerTriangularMatrixCM<double> B({
        {1,0,0},
        {4,5,0},
        {7,8,9},
        {10,11,12}});
    matrix23::LowerTriangularMatrixCM<double> C=A*B;
    EXPECT_EQ(C.row(0),(il{1}));
    EXPECT_EQ(C.row(1),(il{24,30}));
    EXPECT_EQ(C.row(2),(il{77,88,99}));
}
TEST_F(MatrixAlgebraTests, MatrixMultiplyUD)
{
    matrix23::UpperTriangularMatrixCM<double> A({
        {1,2,3,4},
        {0,6,7,8},
        {0,0,11,12}});
    matrix23::DiagonalMatrix<double> B({
        {1,0,0},
        {0,5,0},
        {0,0,9},
        {0,0,0}});
    matrix23::UpperTriangularMatrixCM<double> C=A*B;
    EXPECT_EQ(C.row(0),(il{1,10,27}));
    EXPECT_EQ(C.row(1),(il{30,63}));
    EXPECT_EQ(C.row(2),(il{99}));
  
}
TEST_F(MatrixAlgebraTests, MatrixMultiplyDD)
{
    matrix23::DiagonalMatrix<double> A({
        {1,0,0,0},
        {0,6,0,0},
        {0,0,11,0}});
    matrix23::DiagonalMatrix<double> B({
        {1,0,0},
        {0,5,0},
        {0,0,9},
        {0,0,0}});
    matrix23::DiagonalMatrix<double> C=A*B;
    EXPECT_EQ(C.row(0),(il{1}));
    EXPECT_EQ(C.row(1),(il{30}));
    EXPECT_EQ(C.row(2),(il{99}));
  
}
TEST_F(MatrixAlgebraTests, MatrixMultiplyB1B2)
{
    matrix23::SBandMatrix<double> A({
        {1 , 7, 0, 0, 0, 0},
        {12, 2, 8, 0, 0, 0},
        {0 ,13, 3, 9, 0, 0},
        {0 , 0,14, 4,10, 0},
        {0 , 0, 0,15, 5,11},
        {0 , 0, 0, 0,16, 6}},1);
   
    matrix23::SBandMatrix<double> B({
        {1 , 7,17, 0, 0, 0},
        {12, 2, 8,18, 0, 0},
        {21,13, 3, 9,19, 0},
        {0 ,22,14, 4,10,20},
        {0 , 0,23,15, 5,11},
        {0 , 0, 0,24,16, 6}},2);
    matrix23::SBandMatrix<double> C=A*B;
    EXPECT_EQ(C.packer().bandwidth(),1+2);
    EXPECT_EQ(C.row(0),(il{85,21,73,126}));
    EXPECT_EQ(C.row(1),(il{204,192,244,108,152}));
    EXPECT_EQ(C.row(2),(il{219,263,239,297,147,180}));
    EXPECT_EQ(C.row(3),(il{294,270,328,292,356,190}));
    EXPECT_EQ(C.row(4),(il{330,325,399,351,421}));
    EXPECT_EQ(C.row(5),(il{368,384,176,212}));

}
TEST_F(MatrixAlgebraTests, MatrixMultiplySS)
{
    matrix23::SymmetricMatrixCM<double> A{
        {1,2,3,4},
        {2,5,6,7},
        {3,6,8,9},
        {4,7,9,10}};
    matrix23::SymmetricMatrixCM<double> B{
        {11,12,13,14},
        {12,15,16,17},
        {13,16,18,19},
        {14,17,19,20}};

    matrix23::FullMatrixCM<double> C=A*B;
    EXPECT_EQ(C(0,0),1*11+2*12+3*13+4*14);
    EXPECT_EQ(C(1,0),2*11+5*12+6*13+7*14);
    EXPECT_EQ(C(2,0),3*11+6*12+8*13+9*14);
    EXPECT_EQ(C(3,0),4*11+7*12+9*13+10*14);
    EXPECT_EQ(C(0,1),1*12+2*15+3*16+4*17);
    EXPECT_EQ(C(1,1),2*12+5*15+6*16+7*17);
    EXPECT_EQ(C(2,1),3*12+6*15+8*16+9*17);
    EXPECT_EQ(C(3,1),4*12+7*15+9*16+10*17);
    EXPECT_EQ(C(0,2),1*13+2*16+3*18+4*19);
    EXPECT_EQ(C(1,2),2*13+5*16+6*18+7*19);
    EXPECT_EQ(C(2,2),3*13+6*16+8*18+9*19);
    EXPECT_EQ(C(3,2),4*13+7*16+9*18+10*19);
    EXPECT_EQ(C(0,3),1*14+2*17+3*19+4*20);
    EXPECT_EQ(C(1,3),2*14+5*17+6*19+7*20);
    EXPECT_EQ(C(2,3),3*14+6*17+8*19+9*20);
    EXPECT_EQ(C(3,3),4*14+7*17+9*19+10*20);

    
}
TEST_F(MatrixAlgebraTests, MatrixMultiplySS_commute)
{
    matrix23::SymmetricMatrixCM<double> A{
        {1,2,3,4},
        {2,5,6,7},
        {3,6,8,9},
        {4,7,9,10}};

    matrix23::SymmetricMatrixCM<double> C=A*A; //THis will do a run time check A*A is indeed symmetric.
    EXPECT_EQ(C,(ilil{{30,58,75,85},{58,114,147,167},{75,147,190,216},{85,167,216,246}}));
}
TEST_F(MatrixAlgebraTests, MatrixMultiplySU)
{
    matrix23::SymmetricMatrixCM<double> A{
        {1,2,3,4},
        {2,5,6,7},
        {3,6,8,9},
        {4,7,9,10}};
    matrix23::UpperTriangularMatrixCM<double> B{
        {11,12,13,14},
        { 0,15,16,17},
        { 0, 0,18,19},
        { 0, 0, 0,20}};

    matrix23::FullMatrixCM<double> C=A*B;
    EXPECT_EQ(C,(ilil{{11,42,99,185},{22,99,214,367},{33,126,279,476},{44,153,326,546}}));
    // matrix23::SymmetricMatrixCM<double> D=A*B; Run time fail, A*B is not symmetric.
}
TEST_F(MatrixAlgebraTests, FullColMajor_ops)
{
    matrix23::FullMatrixCM<double> A({
        {1,2,3,4},
        {5,6,7,8},
        {9,10,11,12},
        {13,14,15,16}});

    matrix23::FullMatrixCM<double> C=A*A*A-A+5*A-A*A/2;
    EXPECT_EQ(C,(ilil{{3099,3518,3937,4356},{7187,8142,9097,10052},{11275,12766,14257,15748},{15363,17390,19417,21444}}));
}

TEST_F(MatrixAlgebraTests, ScalarSelfModOps)
{
    matrix23::SymmetricMatrixCM<double> A{
        {1,2,3},
        {2,5,6},
        {3,6,8}};
    matrix23::UpperTriangularMatrixCM<double> B{
        {1,2,3},
        {0,5,6},
        {0,0,8}};
    A+=1;
    EXPECT_EQ(A,(ilil{{2,3,4},{3,6,7},{4,7,9}}));
    B+=1;
    EXPECT_EQ(B,(ilil{{2,3,4},{6,7},{9}}));
    A-=1;
    EXPECT_EQ(A,(ilil{{1,2,3},{2,5,6},{3,6,8}}));
    B-=1;
    EXPECT_EQ(B,(ilil{{1,2,3},{5,6},{8}}));
    A*=-1;
    EXPECT_EQ(A,(ilil{{-1,-2,-3},{-2,-5,-6},{-3,-6,-8}}));
    B*=-1;
    EXPECT_EQ(B,(ilil{{-1,-2,-3},{-5,-6},{-8}}));
    A/=-1;
    EXPECT_EQ(A,(ilil{{1,2,3},{2,5,6},{3,6,8}}));
    B/=-1;
    EXPECT_EQ(B,(ilil{{1,2,3},{5,6},{8}}));
}

TEST_F(MatrixAlgebraTests, MatrixSelfModOps)
{
    matrix23::FullMatrixCM<double> Fcm{3,3,matrix23::one};
    matrix23::FullMatrixRM<double> Frm{3,3,matrix23::one};
    matrix23::SymmetricMatrixCM<double> S{3,matrix23::one};
    matrix23::UpperTriangularMatrixCM<double> U{3,3,matrix23::one};

    Fcm+=Fcm; //use linear iterator
    EXPECT_EQ(Fcm,(ilil{{2,2,2},{2,2,2},{2,2,2}}));
    Fcm-=Fcm;
    EXPECT_EQ(Fcm,(ilil{{0,0,0},{0,0,0},{0,0,0}}));
    Fcm+=Frm; //use op(i,j) or col by col or row by row.
    EXPECT_EQ(Fcm,(ilil{{1,1,1},{1,1,1},{1,1,1}}));
    Fcm-=Frm;
    EXPECT_EQ(Fcm,(ilil{{0,0,0},{0,0,0},{0,0,0}}));
    Fcm+=S; //use op(i,j) or col by col or row by row.
    EXPECT_EQ(Fcm,(ilil{{1,1,1},{1,1,1},{1,1,1}}));
    Fcm-=S;
    EXPECT_EQ(Fcm,(ilil{{0,0,0},{0,0,0},{0,0,0}}));
    Fcm+=U; //use op(i,j) or col by col or row by row.
    EXPECT_EQ(Fcm,(ilil{{1,1,1},{0,1,1},{0,0,1}}));
    Fcm-=U;
    EXPECT_EQ(Fcm,(ilil{{0,0,0},{0,0,0},{0,0,0}}));
    Fcm+=U; //restore

    S+=S;
    EXPECT_EQ(S,(ilil{{2,2,2},{2,2,2},{2,2,2}}));
    S-=S;
    EXPECT_EQ(S,(ilil{{0,0,0},{0,0,0},{0,0,0}}));
    S+=Fcm;
    EXPECT_EQ(S,(ilil{{1,1,1},{1,1,1},{1,1,1}}));
    S-=Fcm;
    EXPECT_EQ(S,(ilil{{0,0,0},{0,0,0},{0,0,0}}));
    S+=U; //use op(i,j) or col by col or row by row.
    EXPECT_EQ(S,(ilil{{1,1,1},{1,1,1},{1,1,1}}));
    S-=U;
    EXPECT_EQ(S,(ilil{{0,0,0},{0,0,0},{0,0,0}}));
    S+=U; //restore

    U+=U;
    EXPECT_EQ(U,(ilil{{2,2,2},{2,2},{2}}));
    U-=U;
    EXPECT_EQ(U,(ilil{{0,0,0},{0,0},{0}}));
    U+=Fcm;
    EXPECT_EQ(U,(ilil{{1,1,1},{1,1},{1}}));
    U-=Fcm;
    EXPECT_EQ(U,(ilil{{0,0,0},{0,0},{0}}));
    U+=S;
    EXPECT_EQ(U,(ilil{{1,1,1},{1,1},{1}}));
    U-=S;
    EXPECT_EQ(U,(ilil{{0,0,0},{0,0},{0}}));

}