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
    v5+=v4;
    EXPECT_EQ(v5,(il{62.5, 257.5, 452.5, 647.5, 842.5}));
    v5-=v4;
    EXPECT_EQ(v5,(il{62, 256.5, 451, 645.5, 840}));
    v5+=1;
    EXPECT_EQ(v5,(il{63, 257.5, 452, 646.5, 841}));
    v5-=1;
    EXPECT_EQ(v5,(il{62, 256.5, 451, 645.5, 840}));
    v5*=-1;
    EXPECT_EQ(v5,(il{-62, -256.5, -451, -645.5, -840}));
    v5/=-1;
    EXPECT_EQ(v5,(il{62, 256.5, 451, 645.5, 840}));

}

#include <valarray>
#include <numeric> // need to include <numeric> for std::ranges::iota! ?
TEST_F(VectorTests, Intrsections)
{
    using matrix23::VectorView;
    using matrix23::intersection;
    using matrix23::inner_product;
    std::valarray<int> v(15);
    std::ranges::iota(v,0); // need to include <numeric> for std::ranges::iota! ?
    auto i1=std::views::iota(size_t(3), size_t(8+1));
    auto i2=std::views::iota(size_t(4), size_t(10+1));
    auto v1v=v | std::views::drop(i1.front()) | std::views::take(i1.size());
    auto v2v=v | std::views::drop(i2.front()) | std::views::take(i2.size());
    auto v1=VectorView(std::move(v1v),i1); // Full view of the valarray
    auto v2=VectorView(std::move(v2v),i2); // Full view of the valarray
  
    intersection inter(i1,i2);
    auto v1_intersection=v1 | std::views::drop(inter.drop1) | std::views::take(inter.indices.size());
    auto v2_intersection=v2 | std::views::drop(inter.drop2) | std::views::take(inter.indices.size());
    int dot=inner_product(v1_intersection,v2_intersection);
    EXPECT_EQ(dot, 4*4 + 5*5 + 6*6 + 7*7 + 8*8); 
}