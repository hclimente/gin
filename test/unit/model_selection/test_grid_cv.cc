//
// Created by hclimente on 30/07/2017.
//

#include "gin/model_selection/grid_cv.h"
#include "gtest/gtest.h"
#include "gin/globals.h"
#include "gin/stats/univariate_association.h"

TEST(GridCV, Constructor) {

	MatrixXd X(10,3);
	X <<    2,0,0,
			2,0,1,
			2,1,0,
			0,1,2,
			1,0,2,
			1,2,0,
			0,2,1,
			0,0,1,
			0,1,0,
			1,0,0;

	VectorXd y(10);
	y <<    0,0,0,0,0,1,1,1,1,2.5;

	MatrixXd dW(4,4);
	dW <<   0, 1, 0, 0,
			1, 0, 1, 0,
			0, 1, 0, 1,
			0, 0, 1, 0;
	SparseMatrixXd W = dW.sparseView();
	VectorXd c(4);
	c << 0.5, 1, 2, 4;

	GridCV g1(&X, &y, &W, 10);
	g1.initFolds(c, SKAT);
	EXPECT_EQ(g1.grids().size(), 10);
	EXPECT_EQ(g1.lambdas().size(), 10);
	EXPECT_NEAR(g1.lambdas().minCoeff(), 0.05, 0.1);
	EXPECT_NEAR(g1.lambdas()(5), 1, 2);
	EXPECT_NEAR(g1.lambdas().maxCoeff(), 40, 0.1);
	EXPECT_NEAR(g1.etas().minCoeff(), 0.5, 0.1);
	EXPECT_NEAR(g1.etas()(5), 1, 2);
	EXPECT_NEAR(g1.etas().maxCoeff(), 4, 0.1);
	EXPECT_EQ(g1.binary_y(), false);

	VectorXd etas(2);
	etas << 4, 8;
	VectorXd lambdas(3);
	lambdas << 3, 5, 7;

	GridCV g2(&X, &y, &W, 10);
	g2.initFolds(etas, lambdas, 0);
	EXPECT_EQ(g2.grids().size(), 10);
	EXPECT_EQ(g2.lambdas().size(), 3);
	EXPECT_NEAR(g2.lambdas().minCoeff(), 3, 0.1);
	EXPECT_NEAR(g2.lambdas().maxCoeff(), 7, 0.1);
	EXPECT_EQ(g2.etas().size(), 2);
	EXPECT_NEAR(g2.etas().minCoeff(), 4, 0.1);
	EXPECT_NEAR(g2.etas().maxCoeff(), 8, 0.1);
	EXPECT_EQ(g2.binary_y(), false);

}

TEST(GridCV, runFolds) {

	MatrixXd X(10, 3);
	X <<    2,0,0,
			2,0,1,
			2,1,0,
			0,1,2,
			1,0,2,
			1,2,0,
			0,2,1,
			0,0,1,
			0,1,0,
			1,0,0;

	VectorXd y(10);
	y <<    0,0,0,0,0,1,1,1,1,1;

	MatrixXd dW(3,3);
	dW <<   0, 1, 0,
			1, 0, 1,
			0, 1, 0;
	SparseMatrixXd W = dW.sparseView();

	VectorXd etas(2);
	etas << 4, 8;
	VectorXd lambdas(3);
	lambdas << 3, 5, 7;

	GridCV g(&X, &y, &W, 10);
	g.initFolds(etas, lambdas, CHI2);

	// k-fold
	// create k-folds
	VectorXd indices = VectorXd::LinSpaced(12,0,12-1);
	Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> perm(12);
	perm.setIdentity();
	indices = perm*indices;

	CCrossValidation cv(0);
	cv.kFold(10, y.rows(), indices);
	g.setCV(cv);

	// values are not initialized
	for (uint i = 0; i < g.grids().size(); i++) {
		for (int e = 0; e < etas.rows(); e++) {
			for (int l = 0; l < lambdas.rows(); l++) {
				double eta = etas[e];
				double lambda = lambdas[l];
				EXPECT_EQ(g.grids()[i] -> selected(eta,lambda).rows(), 0);
			}
		}
	}

	g.runFolds();
	for (uint i = 0; i < g.grids().size(); i++) {
		EXPECT_EQ(g.grids()[0] -> c().rows(), 3);

		for (int e = 0; e < etas.rows(); e++) {
			for (int l = 0; l < lambdas.rows(); l++) {
				double eta = etas[e];
				double lambda = lambdas[l];
				EXPECT_EQ(g.grids()[i] -> selected(eta,lambda).rows(), 3);
			}
		}
	}

	// TODO appropriately split folds, runs them, etc.

}

TEST(GridCV, scoreModels) {

	// only the two first SNPs are causal
	MatrixXd X(12, 54);
	X   <<  2,1,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,
			2,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,
			2,1,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,
			1,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,
			1,2,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,
			1,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2,
			0,0,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2,0,
			0,0,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,
			0,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,
			0,0,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,
			0,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,
			0,0,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2;

	VectorXd y(12);
	y << 1,1,1,1,1,1,0,0,0,0,0,0;

	// network doesnt play a role here
	MatrixXd dW = MatrixXd::Zero(54, 54);
	dW.diagonal(1) = VectorXd::Ones(53);
	dW.diagonal(-1) = VectorXd::Ones(53);
	SparseMatrixXd W = dW.sparseView();
	VectorXd etas(2);
	etas << 0, 4;
	VectorXd lambdas(2);
	lambdas << 0, 1;

	// create k-folds
	VectorXd indices = VectorXd::LinSpaced(12,0,12-1);
	Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> perm(12);
	perm.setIdentity();
	indices = perm*indices;

	CCrossValidation cv(0);
	cv.kFold(10, y.rows(), indices);

	// create and run grids
	GridCV grid_ql(&X, &y, &W, 10);
	grid_ql.initFolds(etas, lambdas, CHI2);
	grid_ql.runFolds();

	grid_ql.scoreModels(AIC);
	EXPECT_EQ(grid_ql.bestEta(), 4);
	EXPECT_NEAR(grid_ql.scoredFolds()(1, 1),
#ifdef __APPLE__
841
#else
821
#endif
	, 1);

	grid_ql.scoreModels(BIC);
	EXPECT_EQ(grid_ql.bestEta(), 4);
	EXPECT_NEAR(grid_ql.scoredFolds()(1, 1),
#ifdef __APPLE__
839
#else
819
#endif
	, 1);

	grid_ql.scoreModels(AICc);
	EXPECT_EQ(grid_ql.bestEta(), 4);
	EXPECT_NEAR(grid_ql.scoredFolds()(1, 1),
#ifdef __APPLE__
838
#else
818
#endif
	, 1);

	grid_ql.scoreModels(CONSISTENCY);
	EXPECT_EQ(grid_ql.bestEta(), 4);
	EXPECT_NEAR(grid_ql.scoredFolds()(1, 1), 1, 0.01);

	// TODO check BIC, AIC, CONSISTENCY for continuous when SKAT is implemented
	// y <<    0.97,0.98,0.94,0.75,0.8,0.77,0.25,0.4,0.25,0.4,0.25,0.3;
	// GridCV grid_qt(&X, &y, &W, etas, lambdas, 10);

}