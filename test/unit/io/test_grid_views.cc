//
// Created by hclimente on 23/08/2017.
//

#include "gin/io/grid_views.h"
#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "mock/model_selection/mock_grid_cv.cc"

using ::testing::Return;
using ::testing::_;
using ::testing::Gt;

TEST(GridViews, ViewSelectedAvg) {

	MockGridCV mockGrid;

	EXPECT_CALL(mockGrid, etas()).WillRepeatedly( Return((VectorXd(3) << 1, 2, 3).finished()) );
	EXPECT_CALL(mockGrid, lambdas()).WillRepeatedly( Return((VectorXd(2) << 4, 5).finished()) );
	EXPECT_CALL(mockGrid, aggregatedFolds(0,0)).WillOnce( Return((VectorXd(10) << 0, 0, 0, 0, 0, 0, 1, 1, 1, 0).finished()) );
	EXPECT_CALL(mockGrid, aggregatedFolds(0,1)).WillOnce( Return((VectorXd(10) << 0, 0, 0, 0, 0, 0, 0, 1, 1, 0).finished()) );
	EXPECT_CALL(mockGrid, aggregatedFolds(1,0)).WillOnce( Return((VectorXd(10) << 0, 0, 0, 0, 1, 1, 1, 1, 1, 1).finished()) );
	EXPECT_CALL(mockGrid, aggregatedFolds(1,1)).WillOnce( Return((VectorXd(10) << 0, 0, 0, 1, 1, 1, 1, 1, 1, 1).finished()) );
	EXPECT_CALL(mockGrid, aggregatedFolds(Gt(1),_)).WillRepeatedly( Return(VectorXd::Zero(10)) );

	GridViews views(&mockGrid);
	MatrixXd expected(4,3);
	expected << 0, 4, 5,
				1, 3, 2,
				2, 6, 7,
				3, 0, 0;
	EXPECT_EQ(views.viewSelectedAvg(), expected);

}

TEST(GridViews, viewSelectionCriterion) {

	MockGridCV mockGrid;

	EXPECT_CALL(mockGrid, etas()).WillRepeatedly( Return((VectorXd(3) << 1, 2, 3).finished()) );
	EXPECT_CALL(mockGrid, lambdas()).WillRepeatedly( Return((VectorXd(2) << 4, 5).finished()) );
	EXPECT_CALL(mockGrid, scoredFolds(0,0)).WillOnce( Return(0.25) );
	EXPECT_CALL(mockGrid, scoredFolds(0,1)).WillOnce( Return(0.6) );
	EXPECT_CALL(mockGrid, scoredFolds(1,0)).WillOnce( Return(0.15) );
	EXPECT_CALL(mockGrid, scoredFolds(1,1)).WillOnce( Return(0.9) );
	EXPECT_CALL(mockGrid, scoredFolds(Gt(1),_)).WillRepeatedly( Return(1) );

	GridViews views(&mockGrid);
	MatrixXd expected(4,3);
	expected << 0, 4, 5,
				1, 0.25, 0.6,
				2, 0.15, 0.9,
				3,    1,   1;
	EXPECT_EQ(views.viewSelectionCriterion(), expected);

}