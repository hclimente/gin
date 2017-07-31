//
// Created by hclimente on 31/07/2017.
//

#include "gtest/gtest.h"
#include "gin/utils/CMatrixHelper.h"
#include "gin/globals.h"

TEST(MatrixHelper, sliceColsMatrixByBinaryVector) {

	MatrixXd X(2,5);
	X << 1,2,3,4,5,
		 6,7,8,9,0;

	MatrixXd sol(2,3);
	sol << 1,3,5,
		   6,8,0;

	VectorXd c(5);
	c << 1, 0, 1, 0, 1;
	EXPECT_EQ(sliceColsMatrixByBinaryVector(X,c), sol);

	c << 0.99, 0.01, 1.01, -0.01, 1;
	EXPECT_EQ(sliceColsMatrixByBinaryVector(X,c), sol);

}
