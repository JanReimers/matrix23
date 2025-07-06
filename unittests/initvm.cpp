// File: initvm.cpp  Test initialization functions for Vector and Matrix.
#include "gtest/gtest.h"
#include <iostream>
#include "matrix23/matrix1.hpp"

using std::cout;
using std::endl;
using matrix23::Vector;
using matrix23::FullMatrixCM;

class InitTests : public ::testing::Test
{
public:
    InitTests() = default;
    ~InitTests() override = default;

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


TEST_F(InitTests, VectorInitializations)
{
    {
        Vector<double> v(10,matrix23::zero);
        EXPECT_EQ(v,(il{0,0,0,0,0,0,0,0,0,0}));
    }
    {
        Vector<double> v(10,matrix23::one);
        EXPECT_EQ(v,(il{1,1,1,1,1,1,1,1,1,1}));
    }
    {
        Vector<double> v(10,matrix23::value);
        EXPECT_EQ(v,(il{1,1,1,1,1,1,1,1,1,1}));
    }
    {
        Vector<double> v(10,matrix23::value,3.14);
        EXPECT_EQ(v,(il{3.14,3.14,3.14,3.14,3.14,3.14,3.14,3.14,3.14,3.14}));
    }
    {
        cout << "These should be random [0,1.0]: ";
        Vector<double> v(10,matrix23::random);
        print(v);
    }   
    {
        cout << "These should be random [0,3.14]: ";
        Vector<double> v(10,matrix23::random,3.14);
        print(v);
    }

    {
        FullMatrixCM<double> m(3,4,matrix23::zero);
        EXPECT_EQ(m,(ilil{{0,0,0,0},{0,0,0,0},{0,0,0,0}}));
    }
    {
        FullMatrixCM<double> m(3,4,matrix23::one);
        EXPECT_EQ(m,(ilil{{1,1,1,1},{1,1,1,1},{1,1,1,1}}));
    }
    {
        FullMatrixCM<double> m(3,4,matrix23::value);
        EXPECT_EQ(m,(ilil{{1,1,1,1},{1,1,1,1},{1,1,1,1}}));
    }
    {
        FullMatrixCM<double> m(3,4,matrix23::value,3.14);
        EXPECT_EQ(m,(ilil{{3.14,3.14,3.14,3.14},{3.14,3.14,3.14,3.14},{3.14,3.14,3.14,3.14}}));
    }
     {
        cout << "These should be random [0,1.0]:";
        FullMatrixCM<double> m(3,4,matrix23::random);
        m.print();
    }   
    {
        cout << "These should be random [0,3.14]:";
        FullMatrixCM<double> m(3,4,matrix23::random,3.14);
        m.print();
    }

}
