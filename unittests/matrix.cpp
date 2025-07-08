// File: matrix.cpp Unit tests for the Matrix<T> class.
#include "gtest/gtest.h"
#include <iostream>
#include <ranges>
#include "matrix23/matrix.hpp"

using std::cout;
using std::endl;
using matrix23::Vector;

namespace matrix23
{
static_assert(isSymmetry<   NoSymmetry<default_data_type<double>,FullPackerCM>> );
static_assert(isSymmetry<    Symmetric<default_data_type<double>,FullPackerCM>> );
static_assert(isSymmetry<AntiSymmetric<default_data_type<double>,FullPackerCM>> );
}

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
    typedef std::initializer_list<il> ilil;
};


TEST_F(MatrixTests, FullColMajor3x4)
{
    matrix23::FullMatrixCM<double> A{
        {1,2,3,4},
        {5,6,7,8},
        {9,10,11,12}};
    EXPECT_EQ(A.row(0),(il{1,2,3,4}));
    EXPECT_EQ(A.row(1),(il{5,6,7,8}));
    EXPECT_EQ(A.row(2),(il{9,10,11,12}));
    EXPECT_EQ(A.col(0),(il{1,5,9}));
    EXPECT_EQ(A.col(1),(il{2,6,10}));
    EXPECT_EQ(A.col(2),(il{3,7,11}));
    EXPECT_EQ(A.col(3),(il{4,8,12}));
}
TEST_F(MatrixTests, FullRowMajor3x4)
{
    matrix23::FullMatrixRM<double> A({
        {1,2,3,4},
        {5,6,7,8},
        {9,10,11,12}});
    EXPECT_EQ(A.row(0),(il{1,2,3,4}));
    EXPECT_EQ(A.row(1),(il{5,6,7,8}));
    EXPECT_EQ(A.row(2),(il{9,10,11,12}));
    EXPECT_EQ(A.col(0),(il{1,5,9}));
    EXPECT_EQ(A.col(1),(il{2,6,10}));
    EXPECT_EQ(A.col(2),(il{3,7,11}));
    EXPECT_EQ(A.col(3),(il{4,8,12}));
}
TEST_F(MatrixTests, FullColMajor4x3)
{
    matrix23::FullMatrixCM<double> A{
        {1,2,3},
        {4,5,6},
        {7,8,9},
        {10,11,12}};
    EXPECT_EQ(A.row(0),(il{1,2,3}));
    EXPECT_EQ(A.row(1),(il{4,5,6}));
    EXPECT_EQ(A.row(2),(il{7,8,9}));
    EXPECT_EQ(A.row(3),(il{10,11,12}));
    EXPECT_EQ(A.col(0),(il{1,4,7,10}));
    EXPECT_EQ(A.col(1),(il{2,5,8,11}));
    EXPECT_EQ(A.col(2),(il{3,6,9,12}));
   
}
TEST_F(MatrixTests, FullRowMajor4x3)
{
    matrix23::FullMatrixRM<double> A({
        {1,2,3},
        {4,5,6},
        {7,8,9},
        {10,11,12}});
    EXPECT_EQ(A.row(0),(il{1,2,3}));
    EXPECT_EQ(A.row(1),(il{4,5,6}));
    EXPECT_EQ(A.row(2),(il{7,8,9}));
    EXPECT_EQ(A.row(3),(il{10,11,12}));
    EXPECT_EQ(A.col(0),(il{1,4,7,10}));
    EXPECT_EQ(A.col(1),(il{2,5,8,11}));
    EXPECT_EQ(A.col(2),(il{3,6,9,12}));
   
}




TEST_F(MatrixTests, UpperTriangularColMajor3x4)
{
    matrix23::UpperTriangularMatrixCM<double> A{
        {1,2,3,4},
        {0,6,7,8},
        {0,0,11,12}};
    EXPECT_EQ(A.row(0),(il{1,2,3,4}));
    EXPECT_EQ(A.row(1),(il{6,7,8}));
    EXPECT_EQ(A.row(2),(il{11,12}));
    EXPECT_EQ(A.col(0),(il{1}));
    EXPECT_EQ(A.col(1),(il{2,6}));
    EXPECT_EQ(A.col(2),(il{3,7,11}));
    EXPECT_EQ(A.col(3),(il{4,8,12}));
}
TEST_F(MatrixTests, UpperTriangularRowMajor3x4)
{
    matrix23::UpperTriangularMatrixRM<double> A({
        {1,2,3,4},
        {0,6,7,8},
        {0,0,11,12}});
    EXPECT_EQ(A.row(0),(il{1,2,3,4}));
    EXPECT_EQ(A.row(1),(il{6,7,8}));
    EXPECT_EQ(A.row(2),(il{11,12}));
    EXPECT_EQ(A.col(0),(il{1}));
    EXPECT_EQ(A.col(1),(il{2,6}));
    EXPECT_EQ(A.col(2),(il{3,7,11}));
    EXPECT_EQ(A.col(3),(il{4,8,12}));
}
TEST_F(MatrixTests, UpperTriangularColMajor4x3)
{
    matrix23::UpperTriangularMatrixCM<double> A{
        {1,2,3},
        {0,5,6},
        {0,0,9},
        {0,0,0}};
    EXPECT_EQ(A.row(0),(il{1,2,3}));
    EXPECT_EQ(A.row(1),(il{5,6}));
    EXPECT_EQ(A.row(2),(il{9}));
    EXPECT_EQ(A.row(3),(il{}));
    EXPECT_EQ(A.col(0),(il{1}));
    EXPECT_EQ(A.col(1),(il{2,5}));
    EXPECT_EQ(A.col(2),(il{3,6,9}));
   
}
TEST_F(MatrixTests, UpperTriangularRowajor4x3)
{
    matrix23::UpperTriangularMatrixRM<double> A({
        {1,2,3},
        {0,5,6},
        {0,0,9},
        {0,0,0}});
    EXPECT_EQ(A.row(0),(il{1,2,3}));
    EXPECT_EQ(A.row(1),(il{5,6}));
    EXPECT_EQ(A.row(2),(il{9}));
    EXPECT_EQ(A.row(3),(il{}));
    EXPECT_EQ(A.col(0),(il{1}));
    EXPECT_EQ(A.col(1),(il{2,5}));
    EXPECT_EQ(A.col(2),(il{3,6,9}));
   
}

TEST_F(MatrixTests, LowerTriangularColMajor3x4)
{
    matrix23::LowerTriangularMatrixCM<double> A{
        {1,0,0,0},
        {5,6,0,0},
        {9,10,11,0}};
    EXPECT_EQ(A.row(0),(il{1}));
    EXPECT_EQ(A.row(1),(il{5,6}));
    EXPECT_EQ(A.row(2),(il{9,10,11}));
    EXPECT_EQ(A.col(0),(il{1,5,9}));
    EXPECT_EQ(A.col(1),(il{6,10}));
    EXPECT_EQ(A.col(2),(il{11}));
    EXPECT_EQ(A.col(3),(il{}));
}
TEST_F(MatrixTests, LowerTriangularRowMajor3x4)
{
    matrix23::LowerTriangularMatrixRM<double> A({
        {1,0,0,0},
        {5,6,0,0},
        {9,10,11,0}});
    EXPECT_EQ(A.row(0),(il{1}));
    EXPECT_EQ(A.row(1),(il{5,6}));
    EXPECT_EQ(A.row(2),(il{9,10,11}));
    EXPECT_EQ(A.col(0),(il{1,5,9}));
    EXPECT_EQ(A.col(1),(il{6,10}));
    EXPECT_EQ(A.col(2),(il{11}));
    EXPECT_EQ(A.col(3),(il{}));
}
TEST_F(MatrixTests, LowerTriangularColMajor4x3)
{
    matrix23::LowerTriangularMatrixCM<double> A{
        {1,0,0},
        {4,5,0},
        {7,8,9},
        {10,11,12}};
    EXPECT_EQ(A.row(0),(il{1}));
    EXPECT_EQ(A.row(1),(il{4,5}));
    EXPECT_EQ(A.row(2),(il{7,8,9}));
    EXPECT_EQ(A.row(3),(il{10,11,12}));
    EXPECT_EQ(A.col(0),(il{1,4,7,10}));
    EXPECT_EQ(A.col(1),(il{5,8,11}));
    EXPECT_EQ(A.col(2),(il{9,12}));
   
}
TEST_F(MatrixTests, LowerTriangularRowMajor4x3)
{
    matrix23::LowerTriangularMatrixRM<double> A({
        {1,0,0},
        {4,5,0},
        {7,8,9},
        {10,11,12}});
    EXPECT_EQ(A.row(0),(il{1}));
    EXPECT_EQ(A.row(1),(il{4,5}));
    EXPECT_EQ(A.row(2),(il{7,8,9}));
    EXPECT_EQ(A.row(3),(il{10,11,12}));
    EXPECT_EQ(A.col(0),(il{1,4,7,10}));
    EXPECT_EQ(A.col(1),(il{5,8,11}));
    EXPECT_EQ(A.col(2),(il{9,12}));
   
}

TEST_F(MatrixTests, Diagonal3x4)
{
    matrix23::DiagonalMatrix<double> A{
        {1,0,0,0},
        {0,6,0,0},
        {0,0,11,0}};
    EXPECT_EQ(A.row(0),(il{1}));
    EXPECT_EQ(A.row(1),(il{6}));
    EXPECT_EQ(A.row(2),(il{11}));
    EXPECT_EQ(A.col(0),(il{1}));
    EXPECT_EQ(A.col(1),(il{6}));
    EXPECT_EQ(A.col(2),(il{11}));
    EXPECT_EQ(A.col(3),(il{}));
}
TEST_F(MatrixTests, Diagonal4x3)
{
    matrix23::DiagonalMatrix<double> A{
        {1,0,0},
        {0,6,0},
        {0,0,11},
        {0,0,0}};
    EXPECT_EQ(A.row(0),(il{1}));
    EXPECT_EQ(A.row(1),(il{6}));
    EXPECT_EQ(A.row(2),(il{11}));
    EXPECT_EQ(A.row(3),(il{}));
    EXPECT_EQ(A.col(0),(il{1}));
    EXPECT_EQ(A.col(1),(il{6}));
    EXPECT_EQ(A.col(2),(il{11}));
}
TEST_F(MatrixTests, SBand6x6_k0)
{
    size_t k=0;
    matrix23::SBandMatrix<double> A({
        {1 , 0, 0, 0, 0, 0},
        {0 , 8, 0, 0, 0, 0},
        {0 , 0,15, 0, 0, 0},
        {0 , 0, 0,22, 0, 0},
        {0 , 0, 0, 0,29, 0},
        {0 , 0, 0, 0, 0,36}
        },k);
    EXPECT_EQ(A.row(0),(il{1}));
    EXPECT_EQ(A.row(1),(il{8}));
    EXPECT_EQ(A.row(2),(il{15}));
    EXPECT_EQ(A.row(3),(il{22}));
    EXPECT_EQ(A.row(4),(il{29}));
    EXPECT_EQ(A.row(5),(il{36}));
    EXPECT_EQ(A.col(0),(il{1}));
    EXPECT_EQ(A.col(1),(il{8}));
    EXPECT_EQ(A.col(2),(il{15}));
    EXPECT_EQ(A.col(3),(il{22}));
    EXPECT_EQ(A.col(4),(il{29}));
    EXPECT_EQ(A.col(5),(il{36}));
}
TEST_F(MatrixTests, SBand6x6_k1)
{
    size_t k=1;
    matrix23::SBandMatrix<double> A({
        {1 ,2 ,0, 0 ,0 ,0},
        {7 ,8 ,9 , 0,0 ,0 },
        {0 ,14,15,16, 0,0 },
        {0 , 0,21,22,23, 0},
        {0 , 0, 0,28,29,30},
        {0 , 0, 0, 0,35,36}
        },k);
    EXPECT_EQ(A.row(0),(il{1,2}));
    EXPECT_EQ(A.row(1),(il{7,8,9}));
    EXPECT_EQ(A.row(2),(il{14,15,16}));
    EXPECT_EQ(A.row(3),(il{21,22,23}));
    EXPECT_EQ(A.row(4),(il{28,29,30}));
    EXPECT_EQ(A.row(5),(il{35,36}));
    EXPECT_EQ(A.col(0),(il{1,7}));
    EXPECT_EQ(A.col(1),(il{2,8,14}));
    EXPECT_EQ(A.col(2),(il{9,15,21}));
    EXPECT_EQ(A.col(3),(il{16,22,28}));
    EXPECT_EQ(A.col(4),(il{23,29,35}));
    EXPECT_EQ(A.col(5),(il{30,36}));
}
TEST_F(MatrixTests, SBand6x6_k2)
{
    size_t k=2;
    matrix23::SBandMatrix<double> A({
        {1 ,2 ,3, 0 ,0 ,0},
        {7 ,8 ,9 ,10,0 ,0 },
        {13,14,15,16,17,0 },
        {0 ,20,21,22,23,24},
        {0 , 0,27,28,29,30},
        {0 , 0, 0,34,35,36}
        },k);
    EXPECT_EQ(A.row(0),(il{1,2,3}));
    EXPECT_EQ(A.row(1),(il{7,8,9,10}));
    EXPECT_EQ(A.row(2),(il{13,14,15,16,17}));
    EXPECT_EQ(A.row(3),(il{20,21,22,23,24}));
    EXPECT_EQ(A.row(4),(il{27,28,29,30}));
    EXPECT_EQ(A.row(5),(il{34,35,36}));
    EXPECT_EQ(A.col(0),(il{1,7,13}));
    EXPECT_EQ(A.col(1),(il{2,8,14,20}));
    EXPECT_EQ(A.col(2),(il{3,9,15,21,27}));
    EXPECT_EQ(A.col(3),(il{10,16,22,28,34}));
    EXPECT_EQ(A.col(4),(il{17,23,29,35}));
    EXPECT_EQ(A.col(5),(il{24,30,36}));
}
TEST_F(MatrixTests, SBand6x6_k3)
{
    size_t k=3;
    matrix23::SBandMatrix<double> A({
        { 1, 2, 3, 4, 0, 0},
        { 7, 8, 9,10,11, 0},
        {13,14,15,16,17,18},
        {19,20,21,22,23,24},
        { 0,26,27,28,29,30},
        { 0, 0,33,34,35,36}
        },k);
    EXPECT_EQ(A.row(0),(il{1,2,3,4}));
    EXPECT_EQ(A.row(1),(il{7,8,9,10,11}));
    EXPECT_EQ(A.row(2),(il{13,14,15,16,17,18}));
    EXPECT_EQ(A.row(3),(il{19,20,21,22,23,24}));
    EXPECT_EQ(A.row(4),(il{26,27,28,29,30}));
    EXPECT_EQ(A.row(5),(il{33,34,35,36}));
    EXPECT_EQ(A.col(0),(il{1,7,13,19}));
    EXPECT_EQ(A.col(1),(il{2,8,14,20,26}));
    EXPECT_EQ(A.col(2),(il{3,9,15,21,27,33}));
    EXPECT_EQ(A.col(3),(il{4,10,16,22,28,34}));
    EXPECT_EQ(A.col(4),(il{11,17,23,29,35}));
    EXPECT_EQ(A.col(5),(il{18,24,30,36}));
}
TEST_F(MatrixTests, SBand6x6_k4)
{
    size_t k=4;
    matrix23::SBandMatrix<double> A({
        { 1, 2, 3, 4, 5, 0},
        { 7, 8, 9,10,11,12},
        {13,14,15,16,17,18},
        {19,20,21,22,23,24},
        {25,26,27,28,29,30},
        { 0,32,33,34,35,36}
        },k);
    EXPECT_EQ(A.row(0),(il{1,2,3,4,5}));
    EXPECT_EQ(A.row(1),(il{7,8,9,10,11,12}));
    EXPECT_EQ(A.row(2),(il{13,14,15,16,17,18}));
    EXPECT_EQ(A.row(3),(il{19,20,21,22,23,24}));
    EXPECT_EQ(A.row(4),(il{25,26,27,28,29,30}));
    EXPECT_EQ(A.row(5),(il{32,33,34,35,36}));
    EXPECT_EQ(A.col(0),(il{1,7,13,19,25}));
    EXPECT_EQ(A.col(1),(il{2,8,14,20,26,32}));
    EXPECT_EQ(A.col(2),(il{3,9,15,21,27,33}));
    EXPECT_EQ(A.col(3),(il{4,10,16,22,28,34}));
    EXPECT_EQ(A.col(4),(il{5,11,17,23,29,35}));
    EXPECT_EQ(A.col(5),(il{12,18,24,30,36}));
}
TEST_F(MatrixTests, SBand6x6_k5)
{
    size_t k=5;
    matrix23::SBandMatrix<double> A({
        { 1, 2, 3, 4, 5, 6},
        { 7, 8, 9,10,11,12},
        {13,14,15,16,17,18},
        {19,20,21,22,23,24},
        {25,26,27,28,29,30},
        {31,32,33,34,35,36}
        },k);
    EXPECT_EQ(A.row(0),(il{1,2,3,4,5,6}));
    EXPECT_EQ(A.row(1),(il{7,8,9,10,11,12}));
    EXPECT_EQ(A.row(2),(il{13,14,15,16,17,18}));
    EXPECT_EQ(A.row(3),(il{19,20,21,22,23,24}));
    EXPECT_EQ(A.row(4),(il{25,26,27,28,29,30}));
    EXPECT_EQ(A.row(5),(il{31,32,33,34,35,36}));
    EXPECT_EQ(A.col(0),(il{1,7,13,19,25,31}));
    EXPECT_EQ(A.col(1),(il{2,8,14,20,26,32}));
    EXPECT_EQ(A.col(2),(il{3,9,15,21,27,33}));
    EXPECT_EQ(A.col(3),(il{4,10,16,22,28,34}));
    EXPECT_EQ(A.col(4),(il{5,11,17,23,29,35}));
    EXPECT_EQ(A.col(5),(il{6,12,18,24,30,36}));
}


TEST_F(MatrixTests, SymmetricColMajor4x4)
{
    matrix23::SymmetricMatrixCM<double> A{
        {1,2,3,4},
        {2,5,6,7},
        {3,6,8,9},
        {4,7,9,10}};
    EXPECT_EQ(A.row(0),(il{1,2,3,4}));
    EXPECT_EQ(A.row(1),(il{2,5,6,7}));
    EXPECT_EQ(A.row(2),(il{3,6,8,9}));
    EXPECT_EQ(A.row(3),(il{4,7,9,10}));
    EXPECT_EQ(A.col(0),(il{1,2,3,4}));
    EXPECT_EQ(A.col(1),(il{2,5,6,7}));
    EXPECT_EQ(A.col(2),(il{3,6,8,9}));
    EXPECT_EQ(A.col(3),(il{4,7,9,10}));
    A(0,2)=0;
    EXPECT_EQ(A.row(0),(il{1,2,0,4}));
    EXPECT_EQ(A.col(2),(il{0,6,8,9}));
}

