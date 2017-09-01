//
// Created by hclimente on 23/07/2017.
//

#include "gin/model_selection/grid.h"
#include "gin/globals.h"

#include "gtest/gtest.h"

// TODO set more setters and getters

TEST(Grid, HasExpectedDimensions) {

	VectorXd c(3);
	c << 1, 1, 1;
	MatrixXd dW(3, 3);
	dW <<   0, 1, 0,
			1, 0, 1,
			0, 1, 0;

	SparseMatrixXd W = dW.sparseView();

	VectorXd etas(4);
	etas << 0, 1, 2, 3;
	VectorXd lambdas(5);
	lambdas << 0, 1, 2, 3, 4;

	std::cout << c << "\n";

	Grid g = Grid(c, &W, etas, lambdas);
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

TEST(Grid, search_runsScones) {

	VectorXd c(10);
	c << 10, 0.4, 0.4, 0.4, 10, 10, 0.4, 0.4, 0.4, 0.4;
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
	VectorXd etas(4);
	etas << 0, 1, 9, 11;
	VectorXd lambdas(2);
	lambdas << 0, 0.3;

	Grid g(c, &W, etas, lambdas);
	g.search();

	VectorXd solution(10);
	solution << 1, 0, 0, 0, 1, 1, 0, 0, 0, 0;
	VectorXd all = VectorXd::Ones(10);
	VectorXd none = VectorXd::Zero(10);

	// same tests run in test_scones
	EXPECT_EQ(g.selected(9,0), solution);
	EXPECT_EQ(g.selected(1,0.3), solution);
	EXPECT_EQ(g.selected(11,0), none);
	EXPECT_EQ(g.selected(0,0), all);

}
