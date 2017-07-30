//
// Created by hclimente on 21/07/2017.
//

#ifndef GIN_UNIVARIATE_ASSOCIATION
#define GIN_UNIVARIATE_ASSOCIATION

#include "gin/stats/CChi2.h"
#include "gin/globals.h"

class UnivariateAssociation {
public:

	// TODO convert to pointers
	UnivariateAssociation(MatrixXd, VectorXd);

	// association scores
	VectorXd computeSKAT();
	VectorXd computeChi2();
	VectorXd computeTrendTest(std::string const&);

private:
	MatrixXd __X;
	VectorXd __y;
	uint __n_features;
	bool __binary_y;
};

#endif //GIN_UNIVARIATE_ASSOCIATION
