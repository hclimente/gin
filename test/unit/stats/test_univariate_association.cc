//
// Created by hclimente on 22/07/2017.
//

#include "gtest/gtest.h"
#include "gin/stats/univariate_association.h"
#include "gin/globals.h"

TEST(testUnivariate, testRightScores2COMPLETE2) {

	MatrixXd tX(2, 60);
	tX <<   0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,
			0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2;

	VectorXd y(60);
	y <<    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;

	UnivariateAssociation u(Eigen::Transpose< MatrixXd >(tX), y);
	VectorXd c = u.computeChi2();

	EXPECT_NEAR(c[0], 10, 0.1);
	EXPECT_NEAR(c[1], 0.4, 0.01);
}

TEST(testUnivariate, testRightScores2COMPLETE) {

	MatrixXd X(10, 5);
	X <<    0, 0, 0, 0, 0,
			0, 0, 1, 1, 1,
			0, 0, 1, 1, 2,
			0, 0, 1, 1, 2,
			0, 0, 1, 1, 2,
			0, 1, 1, 2, 2,
			0, 1, 1, 2, 2,
			0, 1, 1, 2, 2,
			1, 1, 1, 2, 2,
			2, 2, 2, 2, 2;
	VectorXd y(10);
	y << 0, 0, 0, 0, 0, 1, 1, 1, 1, 1;

	UnivariateAssociation u(X, y);

	// chi2
	VectorXd chi2(5);
	chi2 << 2.5, 10, 2, 10, 2.5;

	VectorXd c = u.computeChi2();
	EXPECT_EQ(c.rows(), chi2.rows());

	for (int i = 0; i < 5; i++) {
		EXPECT_NEAR(c[i], chi2[i], 0.1);
	}

	/* trend
	VectorXd trend(5);
	trend << nan, nan, nan, nan, nan;

	c = u.computeTrendTest("additive");
	EXPECT_EQ(c.rows(), trend.rows());

	for (int i = 0; i < 5; i++) {
		EXPECT_NEAR(c[i], trend[i], 0.1);
	}
	*/

	// TODO find right scores
	std::cout << "SKAT\n";

}

