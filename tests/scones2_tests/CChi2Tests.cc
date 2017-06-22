#include "gtest/gtest.h"
#include "CEasyGWAS/stats/CChi2.h"
#include "CEasyGWAS/globals.h"

TEST(CChi2Test, ContingencyTable) {
    VectorXd v1(14);
    VectorXd v2(14);

    v1 << 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5;
    v2 << 2, 2, 2, 1, 2, 2, 1, 2, 1, 2, 2, 1, 2, 2;

    MatrixXd table = CChi2::get2DContingencyTable(v1, v2);

    MatrixXd exp(3,2);
    exp << 1, 4,
           2, 3,
           1, 3;

    ASSERT_EQ(table, exp);

}

TEST(CChi2Test, Chi2Calculation) {

    MatrixXd table(2,3);
    table << 9000, 3000, 1000,
             4500, 10000, 7000;

    float64 chisq = CChi2::calculateChi2(table);
    float64 exp = 8171.02;

    ASSERT_NEAR(chisq, exp, 0.01);

}

TEST(CChi2Test, Chi2TrendCalculation) {

    MatrixXd table(2,3);
    table << 9000, 3000, 1000,
             4500, 10000, 7000;
    VectorXd w(3);
    w << 1, 1, 0;

    float64 T = CChi2::calculateChi2Trend(table, w);
    float64 exp = 2812.4;

    ASSERT_NEAR(T, exp, 0.2);

}

int main(int argc, char* argv[])
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}