//
// Created by hclimente on 22/07/2017.
//

#include "gtest/gtest.h"
#include "gin/stats/univariate_association.h"

TEST(UnivariateAssociation, testComputeChi2) {

	MatrixXd tX(2, 60);
	tX <<   0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,
			0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2;
	MatrixXd X1 = Eigen::Transpose< MatrixXd >(tX);

	VectorXd y1(60);
	y1 <<    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;

	UnivariateAssociation u1(&X1, &y1);
	VectorXd c1 = u1.computeChi2();

	EXPECT_NEAR(c1[0], 9.09, 0.1);
	EXPECT_NEAR(c1[1], 0.36, 0.01);

	MatrixXd X2(10, 5);
	X2 <<   0, 0, 0, 0, 0,
			0, 0, 1, 1, 1,
			0, 0, 1, 1, 2,
			0, 0, 1, 1, 2,
			0, 0, 1, 1, 2,
			0, 1, 1, 2, 2,
			0, 1, 1, 2, 2,
			0, 1, 1, 2, 2,
			1, 1, 1, 2, 2,
			2, 2, 2, 2, 2;
	VectorXd y2(10);
	y2 << 0, 0, 0, 0, 0, 1, 1, 1, 1, 1;

	UnivariateAssociation u2( &X2, &y2);
	VectorXd c2 = u2.computeChi2();

	VectorXd chi2(5);
	chi2 << 1.07, 6.57, 0.67, 6.57, 1.07;

	EXPECT_EQ(c2.rows(), chi2.rows());

	for (int i = 0; i < 5; i++) {
		EXPECT_NEAR(c2[i], chi2[i], 0.01);
	}

}

TEST(UnivariateAssociation, testComputeSKAT) {

	MatrixXd tX(2, 60);
	tX <<   0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,
			0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2;
	MatrixXd X1 = Eigen::Transpose< MatrixXd >(tX);

	VectorXd y1(60);
	y1 <<    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;

	UnivariateAssociation u( &X1, &y1);

	VectorXd s = u.computeSKAT();
	EXPECT_NEAR(s[0], 100, 0.1);
	EXPECT_NEAR(s[1], 4, 0.01);

	MatrixXd X2(10, 5);
	X2 <<    0, 0, 0, 0, 0,
			0, 0, 1, 1, 1,
			0, 0, 1, 1, 2,
			0, 0, 1, 1, 2,
			0, 0, 1, 1, 2,
			0, 1, 1, 2, 2,
			0, 1, 1, 2, 2,
			0, 1, 1, 2, 2,
			1, 1, 1, 2, 2,
			2, 2, 2, 2, 2;
	VectorXd y2(10);
	y2 << 0, 0, 0, 0, 0, 1, 1, 1, 1, 1;

	UnivariateAssociation u2( &X2, &y2);
	VectorXd s2 = u2.computeSKAT();

	VectorXd skat2(5);
	skat2 << 2.25, 9, 1, 9, 2.25;

	EXPECT_EQ(s2.rows(), skat2.rows());

	for (int i = 0; i < 5; i++) {
		EXPECT_NEAR(s2[i], skat2[i], 0.1);
	}


}

