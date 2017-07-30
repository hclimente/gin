//
// Created by hclimente on 21/07/2017.
//

#include "gtest/gtest.h"
#include "gin/feature_selection/scones.h"
#include "gin/globals.h"

TEST(Scones, testHyperparams) {

	VectorXd c(10);
	MatrixXd dW(10, 10);
	dW <<   0, 1, 0, 0, 1, 1, 0, 0, 0, 0,
		    1, 0, 1, 0, 0, 0, 0, 0, 0, 0,
		    0, 1, 0, 1, 0, 0, 0, 0, 0, 0,
		    0, 0, 1, 0, 1, 0, 0, 0, 0, 0,
		    1, 0, 0, 1, 0, 1, 0, 0, 0, 0,
		    1, 0, 0, 0, 1, 0, 1, 0, 0, 0,
		    0, 0, 0, 0, 0, 1, 0, 1, 0, 0,
		    0, 0, 0, 0, 0, 0, 1, 0, 1, 0,
		    0, 0, 0, 0, 0, 0, 0, 1, 0, 1,
		    0, 0, 0, 0, 0, 0, 0, 0, 1, 0;
	SparseMatrixXd W = dW.sparseView();

	VectorXd solution(10);
	solution << 1, 0, 0, 0, 1, 1, 0, 0, 0, 0;

	c << 100, 0, 0, 0, 100, 100, 0, 0, 0, 0;
	Scones etaTest = Scones(c, 1, 0, &W);
	etaTest.selectSnps();

	EXPECT_EQ(etaTest.selected(), solution);
	EXPECT_EQ(etaTest.computeScore(), 297);

	c << 100, 5, 5, 5, 100, 100, 5, 5, 5, 5;
	Scones lambdaTest = Scones(c, 10, 3, &W);
	lambdaTest.selectSnps();

	EXPECT_EQ(lambdaTest.selected(), solution);
	EXPECT_EQ(lambdaTest.computeScore(), 261);

}

TEST(Scones, testHyperparamsTestedInGrids) {

	VectorXd c(10);
	MatrixXd dW(10, 10);
	dW <<   0, 1, 0, 0, 1, 1, 0, 0, 0, 0,
			1, 0, 1, 0, 0, 0, 0, 0, 0, 0,
			0, 1, 0, 1, 0, 0, 0, 0, 0, 0,
			0, 0, 1, 0, 1, 0, 0, 0, 0, 0,
			1, 0, 0, 1, 0, 1, 0, 0, 0, 0,
			1, 0, 0, 0, 1, 0, 1, 0, 0, 0,
			0, 0, 0, 0, 0, 1, 0, 1, 0, 0,
			0, 0, 0, 0, 0, 0, 1, 0, 1, 0,
			0, 0, 0, 0, 0, 0, 0, 1, 0, 1,
			0, 0, 0, 0, 0, 0, 0, 0, 1, 0;
	SparseMatrixXd W = dW.sparseView();

	VectorXd solution(10);
	solution << 1, 0, 0, 0, 1, 1, 0, 0, 0, 0;

	c << 10, 0, 0, 0, 10, 10, 0, 0, 0, 0;
	Scones etaTest = Scones(c, 1, 0, &W);
	etaTest.selectSnps();

	EXPECT_EQ(etaTest.selected(), solution);
	EXPECT_EQ(etaTest.computeScore(), 27);

	c << 10, 0.4, 0.4, 0.4, 10, 10, 0.4, 0.4, 0.4, 0.4;
	Scones lambdaTest = Scones(c, 1, 0.3, &W);
	lambdaTest.selectSnps();

	EXPECT_EQ(lambdaTest.selected(), solution);
	EXPECT_NEAR(lambdaTest.computeScore(), 26.1, 0.01);

}

TEST(Scones, testMaxflow) {

	MatrixXd dW(10, 10);
	dW <<   0, 0, 0, 0, 1, 1, 0, 0, 0, 0,
		    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		    1, 0, 0, 0, 0, 1, 0, 0, 0, 0,
		    1, 0, 0, 0, 1, 0, 0, 0, 0, 0,
		    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		    0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
	SparseMatrixXd W = dW.sparseView();

	VectorXd solution(10);
	solution << 1, 0, 0, 0, 1, 1, 0, 0, 0, 0;
	VectorXd none(10);
	none << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
	VectorXd all(10);
	all << 1, 1, 1, 1, 1, 1, 1, 1, 1, 1;

	// only c != 0 can be connected
	VectorXd c(10);
	c << 100, 0, 0, 0, 100, 100, 0, 0, 0, 0;

	Scones noTerms1 = Scones(c, 0, 0, &W);
	noTerms1.selectSnps();

	EXPECT_EQ(noTerms1.selected(), solution);
	EXPECT_EQ(noTerms1.computeScore(), 300);

	c << 100, 1, 0, 0, 100, 100, 0, 0, 0, 0;

	Scones noTerms2 = Scones(c, 0, 0, &W);
	noTerms2.selectSnps();

	EXPECT_NE(noTerms2.selected(), solution);
	EXPECT_EQ(noTerms2.computeScore(), 301);

	// disconnected from source/sink
	c << 100, 5, 5, 5, 100, 100, 5, 5, 5, 5;
	Scones etaTooLarge = Scones(c, 101, 0, &W);
	etaTooLarge.selectSnps();

	EXPECT_EQ(etaTooLarge.selected(), none);

	Scones etaTooSmall = Scones(c, 3, 0, &W);
	etaTooSmall.selectSnps();

	EXPECT_EQ(etaTooSmall.selected(), all);

}

TEST(Scones, testSettersGetters) {

	VectorXd all(3);
	all << 1, 1, 1;

	VectorXd c(3);
	SparseMatrixXd W;

	Scones empty = Scones(c, 0, 0, &W);
	empty.setSelected(all);

	EXPECT_EQ(empty.selected(), all);

	c << 2, 4, 8;
	Scones empty2 = Scones(c, 0, 0, &W);
	empty2.setSelected(all);

	EXPECT_EQ(empty2.computeScore(), 14);

}

int main(int argc, char* argv[])
{
	testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}