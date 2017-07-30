//
// Created by hclimente on 23/07/2017.
//

#include "gin/model_selection/grid.h"
#include "gtest/gtest.h"
#include "gin/globals.h"

TEST(testGrid, HasExpectedDimensions) {

	MatrixXd X(3, 3);
	X <<    0, 1, 2,
			0, 1, 2,
			0, 1, 2;

	VectorXd y(3);
	y << 0, 0, 1;

	MatrixXd dW(3, 3);
	dW <<   0, 1, 0,
			1, 0, 1,
			0, 1, 0;

	SparseMatrixXd W = dW.sparseView();

	VectorXd etas(4);
	etas << 0, 1, 2, 3;
	VectorXd lambdas(5);
	lambdas << 0, 1, 2, 3, 4;

	Grid g = Grid(X, y, etas, lambdas, &W);
	g.search();

	// check dimensions
	EXPECT_EQ(g.getSelected(etas, lambdas).size(), 20);
	VectorXd oneDimension(1);
	oneDimension << 0;
	EXPECT_EQ(g.getSelected(etas, oneDimension).size(), 4);
	EXPECT_EQ(g.getSelected(oneDimension, lambdas).size(), 5);

	// check all models are run
	std::vector<VectorXd> allModels = g.getSelected(etas, lambdas);

	double max = -1;
	for (VectorXd m : allModels) {
		EXPECT_NE(m.rows(), 0);
	}

}

TEST(testGrid, FindsRightValue) {

	// create dummy scones with preset scores

}

TEST(testGrid, RunsSconesAppropriately) {

	// create dummy scones with preset scores

}
