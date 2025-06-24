// File: vector.cpp Unit tests for the Vector<T> class.
#include "gtest/gtest.h"
#include <iostream>
#include <ranges>
#include "matrix23/vector.hpp"

using std::cout;
using std::endl;
using matrix23::Vector;

class VectorTests : public ::testing::Test
{
public:
    VectorTests() = default;
    ~VectorTests() override = default;

    template <std::ranges::range Range> static void print(Range v)
    {
        cout << "[";
        for (auto element : v) 
            std::cout << element << " ";

        cout << "]\n";
    }
    typedef std::initializer_list<double> il;
};

TEST_F(VectorTests, Initialization)
{
    Vector<double> v1{1, 2, 3, 4, 5};
    EXPECT_EQ(v1.size(), 5);
    EXPECT_EQ(v1,(il{1,2,3,4,5}));
    EXPECT_EQ(v1(0), 1);
    EXPECT_EQ(v1(1), 2);
    EXPECT_EQ(v1(2), 3);
    EXPECT_EQ(v1(3), 4);
    EXPECT_EQ(v1(4), 5);
    // print(v1); // Should print: [1 2 3 4 5 ]
    static_assert(matrix23::isVector<Vector<double>>);
}
TEST_F(VectorTests, Operators)
{
    Vector<double> v1{1, 2, 3, 4, 5};
    // Vector dot product.
    Vector<double> v2{6, 7, 8, 9, 10};
    EXPECT_EQ(v1*v2, 1*6 + 2*7 + 3*8 + 4*9 + 5*10);

    Vector<double> v3 = v1 | std::views::drop(1) | std::views::take(3);
    // print(v3); // [2 3 4 ]
    EXPECT_EQ(v3,(il{2,3,4}));

    Vector<double> v4=v1+v2;
    EXPECT_EQ(v4,(il{7,9,11,13,15}));
    v4=v2-v1; //This works without op=
    EXPECT_EQ(v4,(il{5,5,5,5,5}));
    v4=v1*2;
    EXPECT_EQ(v4,(il{2,4,6,8,10}));
    v4=v1/2;
    EXPECT_EQ(v4,(il{0.5,1,1.5,2,2.5}));

    Vector<double> v5=v4*2+v2-v1*v4*v2+8*v1*v4*v1;
    EXPECT_EQ(v5,(il{62, 256.5, 451, 645.5, 840}));
}