// File: matrix.cpp Unit tests for the Matrix<T> class.
#include "gtest/gtest.h"
#include <iostream>
#include "matrix23/subscriptors.hpp"

using std::cout;
using std::endl;

class PackerTests : public ::testing::Test
{
public:
    PackerTests() = default;
    ~PackerTests() override = default;

    // typedef std::ranges::iota_view<long unsigned int, long unsigned int> iv;
    // friend bool operator==(const iv& a, const iv&b )
    // {
    //     return a.front()==b.front() && a.size()==b.size();
    // }
};

using PackerDeathTest=PackerTests;

using namespace matrix23;

TEST_F(PackerDeathTest, Full)
{
    GTEST_FLAG_SET(death_test_style, "fast"); //Assume single threads in test for now.
    size_t nr=3,nc=4;
    FullPacker fprm(nr,nc,Indexing::row_major);
    EXPECT_EQ(fprm.nr(),3);
    EXPECT_EQ(fprm.nc(),4);
    EXPECT_EQ(fprm.stored_size(),12);
    EXPECT_EQ(fprm.stored_row_size(0),4);
    EXPECT_EQ(fprm.stored_col_size(0),3);
    EXPECT_EQ(fprm.stored_row_size(2),4);
    EXPECT_EQ(fprm.stored_col_size(3),3); 
    EXPECT_EQ(fprm.offset(0,0),0);
    EXPECT_EQ(fprm.offset(1,0),4);
    EXPECT_EQ(fprm.offset(2,0),8);
    EXPECT_EQ(fprm.offset(0,1),1);
    EXPECT_EQ(fprm.offset(1,1),5);
    EXPECT_EQ(fprm.offset(2,3),fprm.stored_size()-1);
#ifdef DO_DEATH_TESTS
    ASSERT_DEATH(fprm.stored_row_size(3),""); //row index 3 is out of bounds.
    ASSERT_DEATH(fprm.stored_col_size(4),""); //col index 4 is out of bounds.
    ASSERT_DEATH(fprm.offset(3,0),"");
    ASSERT_DEATH(fprm.offset(0,4),"");
#endif
    FullPacker fpcm(nr,nc,Indexing::col_major);
    EXPECT_EQ(fpcm.nr(),3);
    EXPECT_EQ(fpcm.nc(),4);
    EXPECT_EQ(fpcm.stored_size(),12);
    EXPECT_EQ(fpcm.stored_row_size(0),4);
    EXPECT_EQ(fpcm.stored_col_size(0),3);
    EXPECT_EQ(fpcm.stored_row_size(2),4);
    EXPECT_EQ(fpcm.stored_col_size(3),3); 
    EXPECT_EQ(fpcm.offset(0,0),0);
    EXPECT_EQ(fpcm.offset(1,0),1);
    EXPECT_EQ(fpcm.offset(2,0),2);
    EXPECT_EQ(fpcm.offset(0,1),3);
    EXPECT_EQ(fpcm.offset(1,1),4);
    EXPECT_EQ(fpcm.offset(2,3),fpcm.stored_size()-1);
#ifdef DO_DEATH_TESTS
    ASSERT_DEATH(fpcm.offset(3,0),"");
    ASSERT_DEATH(fpcm.offset(0,4),"");
    ASSERT_DEATH(fpcm.stored_row_size(3),""); //row index 3 is out of bounds.
    ASSERT_DEATH(fpcm.stored_col_size(4),""); //col index 4 is out of bounds.
#endif
}
TEST_F(PackerDeathTest, UpperTriangular3x4)
{
    GTEST_FLAG_SET(death_test_style, "fast"); //Assume single threads in test for now.
    size_t nr=3,nc=4;
    UpperTriangularPacker utrm(nr,nc,Indexing::row_major);
    EXPECT_EQ(utrm.nr(),3);
    EXPECT_EQ(utrm.nc(),4);
    EXPECT_EQ(utrm.stored_size(),6+3);
    EXPECT_TRUE(utrm.is_stored(0,0));
    EXPECT_TRUE(utrm.is_stored(0,1));
    EXPECT_TRUE(utrm.is_stored(0,2));
    EXPECT_TRUE(utrm.is_stored(0,3));
    EXPECT_TRUE(utrm.is_stored(1,1));
    EXPECT_TRUE(utrm.is_stored(1,2));
    EXPECT_TRUE(utrm.is_stored(1,3));
    EXPECT_TRUE(utrm.is_stored(2,2));
    EXPECT_TRUE(utrm.is_stored(2,3));
    EXPECT_FALSE(utrm.is_stored(1,0));
    EXPECT_FALSE(utrm.is_stored(2,0));
    EXPECT_FALSE(utrm.is_stored(2,1));
    EXPECT_EQ(utrm.stored_row_size(0),4);
    EXPECT_EQ(utrm.stored_row_size(1),3);
    EXPECT_EQ(utrm.stored_row_size(2),2);
    EXPECT_EQ(utrm.stored_col_size(0),1);
    EXPECT_EQ(utrm.stored_col_size(1),2);
    EXPECT_EQ(utrm.stored_col_size(2),3);
    EXPECT_EQ(utrm.stored_col_size(3),3); 
    EXPECT_EQ(utrm.offset(0,0),0);
    EXPECT_EQ(utrm.offset(0,1),1);
    EXPECT_EQ(utrm.offset(1,1),4);
    EXPECT_EQ(utrm.offset(0,2),2);
    EXPECT_EQ(utrm.offset(1,2),5);
    EXPECT_EQ(utrm.offset(2,2),7);
    EXPECT_EQ(utrm.offset(0,3),3);
    EXPECT_EQ(utrm.offset(1,3),6);
    EXPECT_EQ(utrm.offset(2,3),utrm.stored_size()-1);
#ifdef DO_DEATH_TESTS
    ASSERT_DEATH(utrm.stored_row_size(3),""); //row index 3 is out of bounds.
    ASSERT_DEATH(utrm.stored_col_size(4),""); //col index 4 is out of bounds.
    ASSERT_DEATH(utrm.offset(3,0),"");
    ASSERT_DEATH(utrm.offset(0,4),"");
#endif
    UpperTriangularPacker utcm(nr,nc,Indexing::col_major);
    EXPECT_EQ(utcm.nr(),3);
    EXPECT_EQ(utcm.nc(),4);
    EXPECT_EQ(utcm.stored_size(),6+3);
    EXPECT_TRUE(utcm.is_stored(0,0));
    EXPECT_TRUE(utcm.is_stored(0,1));
    EXPECT_TRUE(utcm.is_stored(0,2));
    EXPECT_TRUE(utcm.is_stored(0,3));
    EXPECT_TRUE(utcm.is_stored(1,1));
    EXPECT_TRUE(utcm.is_stored(1,2));
    EXPECT_TRUE(utcm.is_stored(1,3));
    EXPECT_TRUE(utcm.is_stored(2,2));
    EXPECT_TRUE(utcm.is_stored(2,3));
    EXPECT_FALSE(utcm.is_stored(1,0));
    EXPECT_FALSE(utcm.is_stored(2,0));
    EXPECT_FALSE(utcm.is_stored(2,1));
    EXPECT_EQ(utcm.stored_row_size(0),4);
    EXPECT_EQ(utcm.stored_row_size(1),3);
    EXPECT_EQ(utcm.stored_row_size(2),2);
    EXPECT_EQ(utcm.stored_col_size(0),1);
    EXPECT_EQ(utcm.stored_col_size(1),2);
    EXPECT_EQ(utcm.stored_col_size(2),3);
    EXPECT_EQ(utcm.stored_col_size(3),3); 
    EXPECT_EQ(utcm.offset(0,0),0);
    EXPECT_EQ(utcm.offset(0,1),1);
    EXPECT_EQ(utcm.offset(1,1),2);
    EXPECT_EQ(utcm.offset(0,2),3);
    EXPECT_EQ(utcm.offset(1,2),4); 
    EXPECT_EQ(utcm.offset(2,2),5);
    EXPECT_EQ(utcm.offset(0,3),6);
    EXPECT_EQ(utcm.offset(1,3),7);
    EXPECT_EQ(utcm.offset(2,3),utcm.stored_size()-1);
#ifdef DO_DEATH_TESTS
    ASSERT_DEATH(utcm.offset(3,0),"");
    ASSERT_DEATH(utcm.offset(0,4),"");
    ASSERT_DEATH(utcm.stored_row_size(3),""); //row index 3 is out of bounds.
    ASSERT_DEATH(utcm.stored_col_size(4),""); //col index 4 is out of bounds.
#endif
}
TEST_F(PackerDeathTest, UpperTriangular4x3)
{
    GTEST_FLAG_SET(death_test_style, "fast"); //Assume single threads in test for now.
    size_t nr=4,nc=3;
    UpperTriangularPacker utrm(nr,nc,Indexing::row_major);
    EXPECT_EQ(utrm.nr(),4);
    EXPECT_EQ(utrm.nc(),3);
    EXPECT_EQ(utrm.stored_size(),6);
    EXPECT_TRUE(utrm.is_stored(0,0));
    EXPECT_TRUE(utrm.is_stored(0,1));
    EXPECT_TRUE(utrm.is_stored(0,2));
    EXPECT_TRUE(utrm.is_stored(1,1));
    EXPECT_TRUE(utrm.is_stored(1,2));
    EXPECT_TRUE(utrm.is_stored(2,2));
    EXPECT_FALSE(utrm.is_stored(1,0));
    EXPECT_FALSE(utrm.is_stored(2,0));
    EXPECT_FALSE(utrm.is_stored(2,1));
    EXPECT_FALSE(utrm.is_stored(3,0));
    EXPECT_FALSE(utrm.is_stored(3,1));
    EXPECT_FALSE(utrm.is_stored(3,2));
    EXPECT_EQ(utrm.stored_row_size(0),3);
    EXPECT_EQ(utrm.stored_row_size(1),2);
    EXPECT_EQ(utrm.stored_row_size(2),1);
    EXPECT_EQ(utrm.stored_col_size(0),1);
    EXPECT_EQ(utrm.stored_col_size(1),2);
    EXPECT_EQ(utrm.stored_col_size(2),3);
    EXPECT_EQ(utrm.offset(0,0),0);
    EXPECT_EQ(utrm.offset(0,1),1);
    EXPECT_EQ(utrm.offset(1,1),3);
    EXPECT_EQ(utrm.offset(0,2),2);
    EXPECT_EQ(utrm.offset(1,2),4);
    EXPECT_EQ(utrm.offset(2,2),5);

    UpperTriangularPacker utcm(nr,nc,Indexing::col_major);
    EXPECT_EQ(utcm.nr(),4);
    EXPECT_EQ(utcm.nc(),3);
    EXPECT_EQ(utcm.stored_size(),6);
    EXPECT_TRUE(utcm.is_stored(0,0));
    EXPECT_TRUE(utcm.is_stored(0,1));
    EXPECT_TRUE(utcm.is_stored(0,2));
    EXPECT_TRUE(utcm.is_stored(1,1));
    EXPECT_TRUE(utcm.is_stored(1,2));
    EXPECT_TRUE(utcm.is_stored(2,2));
    EXPECT_FALSE(utcm.is_stored(1,0));
    EXPECT_FALSE(utcm.is_stored(2,0));
    EXPECT_FALSE(utcm.is_stored(2,1));
    EXPECT_FALSE(utcm.is_stored(3,0));
    EXPECT_FALSE(utcm.is_stored(3,1));
    EXPECT_FALSE(utcm.is_stored(3,2));

    EXPECT_EQ(utcm.stored_row_size(0),3);
    EXPECT_EQ(utcm.stored_row_size(1),2);
    EXPECT_EQ(utcm.stored_row_size(2),1);
    EXPECT_EQ(utcm.stored_row_size(3),0); 
    EXPECT_EQ(utcm.stored_col_size(0),1);
    EXPECT_EQ(utcm.stored_col_size(1),2);
    EXPECT_EQ(utcm.stored_col_size(2),3);
    EXPECT_EQ(utcm.offset(0,0),0);
    EXPECT_EQ(utcm.offset(0,1),1);
    EXPECT_EQ(utcm.offset(1,1),2);
    EXPECT_EQ(utcm.offset(0,2),3);
    EXPECT_EQ(utcm.offset(1,2),4); 
    EXPECT_EQ(utcm.offset(2,2),5);
}
TEST_F(PackerDeathTest, LowerTriangular3x4)
{
    GTEST_FLAG_SET(death_test_style, "fast"); //Assume single threads in test for now.
    size_t nr=3,nc=4;
    LowerTriangularPacker ltrm(nr,nc,Indexing::row_major);
    EXPECT_EQ(ltrm.nr(),3);
    EXPECT_EQ(ltrm.nc(),4);
    EXPECT_EQ(ltrm.stored_size(),6);
    EXPECT_TRUE(ltrm.is_stored(0,0));
    EXPECT_TRUE(ltrm.is_stored(1,0));
    EXPECT_TRUE(ltrm.is_stored(2,0));
    EXPECT_TRUE(ltrm.is_stored(1,1));
    EXPECT_TRUE(ltrm.is_stored(2,1));
    EXPECT_TRUE(ltrm.is_stored(2,2));
    EXPECT_FALSE(ltrm.is_stored(0,1));
    EXPECT_FALSE(ltrm.is_stored(0,2));
    EXPECT_FALSE(ltrm.is_stored(1,2));
    EXPECT_FALSE(ltrm.is_stored(0,3));
    EXPECT_FALSE(ltrm.is_stored(1,3));
    EXPECT_FALSE(ltrm.is_stored(2,3));

    EXPECT_EQ(ltrm.stored_row_size(0),1);
    EXPECT_EQ(ltrm.stored_row_size(1),2);
    EXPECT_EQ(ltrm.stored_row_size(2),3);
    EXPECT_EQ(ltrm.stored_col_size(0),3);
    EXPECT_EQ(ltrm.stored_col_size(1),2);
    EXPECT_EQ(ltrm.stored_col_size(2),1);
    EXPECT_EQ(ltrm.stored_col_size(3),0); 
    EXPECT_EQ(ltrm.offset(0,0),0);
    EXPECT_EQ(ltrm.offset(1,0),1);
    EXPECT_EQ(ltrm.offset(1,1),2);
    EXPECT_EQ(ltrm.offset(2,0),3);
    EXPECT_EQ(ltrm.offset(2,1),4);
    EXPECT_EQ(ltrm.offset(2,2),5);

    LowerTriangularPacker ltcm(nr,nc,Indexing::col_major);
    EXPECT_EQ(ltcm.nr(),3);
    EXPECT_EQ(ltcm.nc(),4);
    EXPECT_EQ(ltcm.stored_size(),6);
    EXPECT_TRUE(ltcm.is_stored(0,0));
    EXPECT_TRUE(ltcm.is_stored(1,0));
    EXPECT_TRUE(ltcm.is_stored(2,0));
    EXPECT_TRUE(ltcm.is_stored(1,1));
    EXPECT_TRUE(ltcm.is_stored(2,1));
    EXPECT_TRUE(ltcm.is_stored(2,2));
    EXPECT_FALSE(ltcm.is_stored(0,1));
    EXPECT_FALSE(ltcm.is_stored(0,2));
    EXPECT_FALSE(ltcm.is_stored(1,2));
    EXPECT_FALSE(ltcm.is_stored(0,3));
    EXPECT_FALSE(ltcm.is_stored(1,3));
    EXPECT_FALSE(ltcm.is_stored(2,3));

    EXPECT_EQ(ltcm.stored_row_size(0),1);
    EXPECT_EQ(ltcm.stored_row_size(1),2);
    EXPECT_EQ(ltcm.stored_row_size(2),3);
    EXPECT_EQ(ltcm.stored_col_size(0),3);
    EXPECT_EQ(ltcm.stored_col_size(1),2);
    EXPECT_EQ(ltcm.stored_col_size(2),1);
    EXPECT_EQ(ltcm.stored_col_size(3),0); 
    EXPECT_EQ(ltcm.offset(0,0),0);
    EXPECT_EQ(ltcm.offset(1,0),1);
    EXPECT_EQ(ltcm.offset(1,1),3);
    EXPECT_EQ(ltcm.offset(2,0),2);
    EXPECT_EQ(ltcm.offset(2,1),4); 
    EXPECT_EQ(ltcm.offset(2,2),5); 

}
TEST_F(PackerDeathTest, LowerTriangular4x3)
{
    GTEST_FLAG_SET(death_test_style, "fast"); //Assume single threads in test for now.
    size_t nr=4,nc=3;
    LowerTriangularPacker ltrm(nr,nc,Indexing::row_major);
    EXPECT_EQ(ltrm.nr(),4);
    EXPECT_EQ(ltrm.nc(),3);
    EXPECT_EQ(ltrm.stored_size(),6+3);
    EXPECT_TRUE(ltrm.is_stored(0,0));
    EXPECT_TRUE(ltrm.is_stored(1,0));
    EXPECT_TRUE(ltrm.is_stored(2,0));
    EXPECT_TRUE(ltrm.is_stored(1,1));
    EXPECT_TRUE(ltrm.is_stored(2,1));
    EXPECT_TRUE(ltrm.is_stored(2,2));
    EXPECT_TRUE(ltrm.is_stored(3,0));
    EXPECT_TRUE(ltrm.is_stored(3,1));
    EXPECT_TRUE(ltrm.is_stored(3,2));
    EXPECT_FALSE(ltrm.is_stored(0,1));
    EXPECT_FALSE(ltrm.is_stored(0,2));
    EXPECT_FALSE(ltrm.is_stored(1,2));

    EXPECT_EQ(ltrm.stored_row_size(0),1);
    EXPECT_EQ(ltrm.stored_row_size(1),2);
    EXPECT_EQ(ltrm.stored_row_size(2),3);
    EXPECT_EQ(ltrm.stored_row_size(3),3); 
    EXPECT_EQ(ltrm.stored_col_size(0),4);
    EXPECT_EQ(ltrm.stored_col_size(1),3);
    EXPECT_EQ(ltrm.stored_col_size(2),2);
    EXPECT_EQ(ltrm.offset(0,0),0);
    EXPECT_EQ(ltrm.offset(1,0),1);
    EXPECT_EQ(ltrm.offset(1,1),2);
    EXPECT_EQ(ltrm.offset(2,0),3);
    EXPECT_EQ(ltrm.offset(2,1),4);
    EXPECT_EQ(ltrm.offset(2,2),5);
    EXPECT_EQ(ltrm.offset(3,0),6);
    EXPECT_EQ(ltrm.offset(3,1),7);
    EXPECT_EQ(ltrm.offset(3,2),8);

    LowerTriangularPacker ltcm(nr,nc,Indexing::col_major);
    EXPECT_EQ(ltcm.nr(),4);
    EXPECT_EQ(ltcm.nc(),3);
    EXPECT_EQ(ltcm.stored_size(),6+3);
    EXPECT_TRUE(ltcm.is_stored(0,0));
    EXPECT_TRUE(ltcm.is_stored(1,0));
    EXPECT_TRUE(ltcm.is_stored(2,0));
    EXPECT_TRUE(ltcm.is_stored(1,1));
    EXPECT_TRUE(ltcm.is_stored(2,1));
    EXPECT_TRUE(ltcm.is_stored(2,2));
    EXPECT_TRUE(ltcm.is_stored(3,0));
    EXPECT_TRUE(ltcm.is_stored(3,1));
    EXPECT_TRUE(ltcm.is_stored(3,2));
    EXPECT_FALSE(ltcm.is_stored(0,1));
    EXPECT_FALSE(ltcm.is_stored(0,2));
    EXPECT_FALSE(ltcm.is_stored(1,2));
    

    EXPECT_EQ(ltcm.stored_row_size(0),1);
    EXPECT_EQ(ltcm.stored_row_size(1),2);
    EXPECT_EQ(ltcm.stored_row_size(2),3);
    EXPECT_EQ(ltcm.stored_row_size(3),3);
    EXPECT_EQ(ltcm.stored_col_size(0),4);
    EXPECT_EQ(ltcm.stored_col_size(1),3);
    EXPECT_EQ(ltcm.stored_col_size(2),2);
    EXPECT_EQ(ltcm.offset(0,0),0);
    EXPECT_EQ(ltcm.offset(1,0),1);
    EXPECT_EQ(ltcm.offset(2,0),2);
    EXPECT_EQ(ltcm.offset(3,0),3);
    EXPECT_EQ(ltcm.offset(1,1),4);
    EXPECT_EQ(ltcm.offset(2,1),5); 
    EXPECT_EQ(ltcm.offset(3,1),6); 
    EXPECT_EQ(ltcm.offset(2,2),7); 
    EXPECT_EQ(ltcm.offset(3,2),8); 

}

TEST_F(PackerDeathTest, DiagonalPacker3x4)
{
    size_t nr=3,nc=4;
    DiagonalPacker d(nr,nc);
    EXPECT_EQ(d.nr(),3);
    EXPECT_EQ(d.nc(),4);
    EXPECT_EQ(d.stored_size(),3);
    EXPECT_EQ(d.stored_row_size(0),1);
    EXPECT_EQ(d.stored_row_size(1),1);
    EXPECT_EQ(d.stored_row_size(2),1);
    EXPECT_EQ(d.stored_col_size(0),1);
    EXPECT_EQ(d.stored_col_size(1),1);
    EXPECT_EQ(d.stored_col_size(2),1);
    EXPECT_EQ(d.stored_col_size(3),0);
    EXPECT_TRUE(d.is_stored(0,0)); 
    EXPECT_TRUE(d.is_stored(1,1)); 
    EXPECT_TRUE(d.is_stored(2,2)); 
    EXPECT_FALSE(d.is_stored(1,0)); 
    EXPECT_FALSE(d.is_stored(0,1)); 
    EXPECT_FALSE(d.is_stored(0,2)); 
    EXPECT_FALSE(d.is_stored(2,0)); 
    EXPECT_EQ(d.offset(0,0),0);
    EXPECT_EQ(d.offset(1,1),1);
    EXPECT_EQ(d.offset(2,2),2);
}

TEST_F(PackerDeathTest, SBandPacker6x6)
{
    size_t n=6,k=2;
    SBandPacker sb(n,k);
    EXPECT_EQ(sb.nr(),6);
    EXPECT_EQ(sb.nc(),6);
    EXPECT_EQ(sb.stored_size(),6*5);
    EXPECT_EQ(sb.stored_row_size(0),3);
    EXPECT_EQ(sb.stored_row_size(1),4);
    EXPECT_EQ(sb.stored_row_size(2),5);
    EXPECT_EQ(sb.stored_row_size(3),5);
    EXPECT_EQ(sb.stored_row_size(4),4);
    EXPECT_EQ(sb.stored_row_size(5),3);
    EXPECT_EQ(sb.stored_col_size(0),3);
    EXPECT_EQ(sb.stored_col_size(1),4);
    EXPECT_EQ(sb.stored_col_size(2),5);
    EXPECT_EQ(sb.stored_col_size(3),5);
    EXPECT_EQ(sb.stored_col_size(4),4);
    EXPECT_EQ(sb.stored_col_size(5),3);
    EXPECT_TRUE(sb.is_stored(0,0)); 
    EXPECT_TRUE(sb.is_stored(0,1)); 
    EXPECT_TRUE(sb.is_stored(0,2)); 
    EXPECT_TRUE(sb.is_stored(1,0)); 
    EXPECT_TRUE(sb.is_stored(1,1)); 
    EXPECT_TRUE(sb.is_stored(1,2)); 
    EXPECT_TRUE(sb.is_stored(1,3)); 
    EXPECT_TRUE(sb.is_stored(2,0)); 
    EXPECT_TRUE(sb.is_stored(2,1)); 
    EXPECT_TRUE(sb.is_stored(2,2)); 
    EXPECT_TRUE(sb.is_stored(2,3)); 
    EXPECT_TRUE(sb.is_stored(2,4)); 
    EXPECT_TRUE(sb.is_stored(3,1)); 
    EXPECT_TRUE(sb.is_stored(3,2)); 
    EXPECT_TRUE(sb.is_stored(3,3)); 
    EXPECT_TRUE(sb.is_stored(3,4)); 
    EXPECT_TRUE(sb.is_stored(3,5)); 
    EXPECT_TRUE(sb.is_stored(4,2)); 
    EXPECT_TRUE(sb.is_stored(4,3)); 
    EXPECT_TRUE(sb.is_stored(4,4)); 
    EXPECT_TRUE(sb.is_stored(4,5)); 
    EXPECT_TRUE(sb.is_stored(5,3)); 
    EXPECT_TRUE(sb.is_stored(5,4)); 
    EXPECT_TRUE(sb.is_stored(5,5)); 
    EXPECT_FALSE(sb.is_stored(0,3)); 
    EXPECT_FALSE(sb.is_stored(0,4)); 
    EXPECT_FALSE(sb.is_stored(0,5)); 
    EXPECT_FALSE(sb.is_stored(1,4)); 
    EXPECT_FALSE(sb.is_stored(1,5)); 
    EXPECT_FALSE(sb.is_stored(2,5)); 
    EXPECT_FALSE(sb.is_stored(3,0)); 
    EXPECT_FALSE(sb.is_stored(4,0)); 
    EXPECT_FALSE(sb.is_stored(4,1)); 
    EXPECT_FALSE(sb.is_stored(5,0)); 
    EXPECT_FALSE(sb.is_stored(5,1)); 
    EXPECT_FALSE(sb.is_stored(5,2)); 
    EXPECT_EQ(sb.offset(0,0),2);
    EXPECT_EQ(sb.offset(1,0),3);
    EXPECT_EQ(sb.offset(2,0),4);
    EXPECT_EQ(sb.offset(0,1),6);
    EXPECT_EQ(sb.offset(1,1),7);
    EXPECT_EQ(sb.offset(2,1),8);
    EXPECT_EQ(sb.offset(3,1),9 );
    EXPECT_EQ(sb.offset(0,2),10);
    EXPECT_EQ(sb.offset(1,2),11);
    EXPECT_EQ(sb.offset(2,2),12);
    EXPECT_EQ(sb.offset(3,2),13);
    EXPECT_EQ(sb.offset(4,2),14);
    EXPECT_EQ(sb.offset(1,3),15);
    EXPECT_EQ(sb.offset(2,3),16);
    EXPECT_EQ(sb.offset(3,3),17);
    EXPECT_EQ(sb.offset(4,3),18);
    EXPECT_EQ(sb.offset(5,3),19);
    EXPECT_EQ(sb.offset(2,4),20);
    EXPECT_EQ(sb.offset(3,4),21);
    EXPECT_EQ(sb.offset(4,4),22);
    EXPECT_EQ(sb.offset(5,4),23);
    EXPECT_EQ(sb.offset(3,5),25);
    EXPECT_EQ(sb.offset(4,5),26);
    EXPECT_EQ(sb.offset(5,5),27);
    
}

namespace std::ranges {
 
    bool operator==(const iota_view<size_t,size_t>& a, const iota_view<size_t,size_t>&b )
    {
        return a.front()==b.front() && a.size()==b.size();
    }
}

TEST_F(PackerTests,FullShaper3x4)
{
    size_t nr=3, nc=4;
    FullShaper fs(nr,nc);
    for (size_t col=0;col<nc;col++)
        EXPECT_EQ(fs.nonzero_row_indexes(col),iota_view(0,3));
    for (size_t row=0;row<nr;row++)
        EXPECT_EQ(fs.nonzero_col_indexes(row),iota_view(0,4));
}
TEST_F(PackerTests,FullShaper4x3)
{
    size_t nr=4, nc=3;
    FullShaper fs(nr,nc);
    for (size_t col=0;col<nc;col++)
        EXPECT_EQ(fs.nonzero_row_indexes(col),iota_view(0,4));
    for (size_t row=0;row<nr;row++)
        EXPECT_EQ(fs.nonzero_col_indexes(row),iota_view(0,3));
}
TEST_F(PackerTests,UpperTriangularShaper3x4)
{
    size_t nr=3, nc=4;
    UpperTriangularShaper s(nr,nc);

    EXPECT_EQ(s.nonzero_row_indexes(0),iota_view(0,1));
    EXPECT_EQ(s.nonzero_row_indexes(1),iota_view(0,2));
    EXPECT_EQ(s.nonzero_row_indexes(2),iota_view(0,3));
    EXPECT_EQ(s.nonzero_row_indexes(3),iota_view(0,3));
    EXPECT_EQ(s.nonzero_col_indexes(0),iota_view(0,4));
    EXPECT_EQ(s.nonzero_col_indexes(1),iota_view(1,4));
    EXPECT_EQ(s.nonzero_col_indexes(2),iota_view(2,4));
}

TEST_F(PackerTests,UpperTriangularShaper4x3)
{
    size_t nr=4, nc=3;
    UpperTriangularShaper s(nr,nc);

    EXPECT_EQ(s.nonzero_row_indexes(0),iota_view(0,1));
    EXPECT_EQ(s.nonzero_row_indexes(1),iota_view(0,2));
    EXPECT_EQ(s.nonzero_row_indexes(2),iota_view(0,3));
    EXPECT_EQ(s.nonzero_col_indexes(0),iota_view(0,3));
    EXPECT_EQ(s.nonzero_col_indexes(1),iota_view(1,3));
    EXPECT_EQ(s.nonzero_col_indexes(2),iota_view(2,3));
    EXPECT_TRUE(s.nonzero_col_indexes(3).empty());
}
TEST_F(PackerTests,LowerTriangularShaper3x4)
{
    size_t nr=3, nc=4;
    LowerTriangularShaper s(nr,nc);

    EXPECT_EQ(s.nonzero_row_indexes(0),iota_view(0,3));
    EXPECT_EQ(s.nonzero_row_indexes(1),iota_view(1,3));
    EXPECT_EQ(s.nonzero_row_indexes(2),iota_view(2,3));
    EXPECT_TRUE(s.nonzero_row_indexes(3).empty());
    EXPECT_EQ(s.nonzero_col_indexes(0),iota_view(0,1));
    EXPECT_EQ(s.nonzero_col_indexes(1),iota_view(0,2));
    EXPECT_EQ(s.nonzero_col_indexes(2),iota_view(0,3));
}
TEST_F(PackerTests,LowerTriangularShaper4x3)
{
    size_t nr=4, nc=3;
    LowerTriangularShaper s(nr,nc);

    EXPECT_EQ(s.nonzero_row_indexes(0),iota_view(0,4));
    EXPECT_EQ(s.nonzero_row_indexes(1),iota_view(1,4));
    EXPECT_EQ(s.nonzero_row_indexes(2),iota_view(2,4));
    EXPECT_EQ(s.nonzero_col_indexes(0),iota_view(0,1));
    EXPECT_EQ(s.nonzero_col_indexes(1),iota_view(0,2));
    EXPECT_EQ(s.nonzero_col_indexes(2),iota_view(0,3));
    EXPECT_EQ(s.nonzero_col_indexes(3),iota_view(0,3));

}
