// File: unittests/blas.cpp  Run some tests using the blas library.

#include "gtest/gtest.h"
#include <iostream>
#include "matrix23/matrix1.hpp"
#include "matrix23/blas.hpp"

using std::cout;
using std::endl;
using matrix23::Vector;
using matrix23::FullMatrixCM;

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
    FullMatrixCM<double> A(nr,nc,matrix23::random);
    Vector<double> x(nc,matrix23::random),y(nr);
    matrix23::gemv(1.0,A,x,0.0,y);
    EXPECT_EQ(A*x,y);

}