#include "unit_test.h"

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    ::testing::GTEST_FLAG(filter) = "GlobalTest.min_max_3_variable_quadratic_polynomial";
    int result = RUN_ALL_TESTS();
    return result;
}