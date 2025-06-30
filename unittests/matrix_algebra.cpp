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
