//
// Created by hclimente on 23/07/2017.
//

#include "gin/model_selection/grid.h"
#include "gtest/gtest.h"
#include "gin/globals.h"

TEST(Grid, HasExpectedDimensions) {

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

	Grid g = Grid(X, y, &W, CHI2, etas, lambdas);
	g.search();

	// check dimensions
	EXPECT_EQ(g.selected(etas, lambdas).size(), 20);
	VectorXd oneDimension(1);
	oneDimension << 0;
	EXPECT_EQ(g.selected(etas, oneDimension).size(), 4);
	EXPECT_EQ(g.selected(oneDimension, lambdas).size(), 5);

	// check all models are run
	std::vector<VectorXd> allModels = g.selected(etas, lambdas);

	for (VectorXd m : allModels) {
		EXPECT_NE(m.rows(), 0);
	}

}

TEST(Grid, RunsScones) {

	// transposed matrix
	// chi2 equals to 10 on associated snps and to 0.4 on unassociated snps
	MatrixXd tX(10, 60);
	tX <<   0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,
			0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,
			0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,
			0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,
			0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,
			0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,
			0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,
			0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,
			0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,
			0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2;

	VectorXd y(60);
	y <<    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;

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
	VectorXd all(10);
	all << 1, 1, 1, 1, 1, 1, 1, 1, 1, 1;
	VectorXd none(10);
	none << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;

	VectorXd etas(4);
	etas << 0, 1, 10, 11;
	VectorXd lambdas(2);
	lambdas << 0, 0.3;

	Grid g(Eigen::Transpose< MatrixXd >(tX), y, &W, 1, etas, lambdas);
	g.search();

	// same tests run in test_scones
	EXPECT_EQ(g.selected(10,0), solution);
	EXPECT_EQ(g.selected(1,0.3), solution);
	EXPECT_EQ(g.selected(11,0), none);
	EXPECT_EQ(g.selected(0,0), all);

}
