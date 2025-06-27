#include "gtest/gtest.h"

int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    GTEST_FLAG_SET(death_test_style, "fast"); //Assume single threads in test for now.
    return RUN_ALL_TESTS();
}


