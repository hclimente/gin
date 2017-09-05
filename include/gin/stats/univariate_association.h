//
// Created by hclimente on 21/07/2017.
//

#ifndef GIN_UNIVARIATE_ASSOCIATION
#define GIN_UNIVARIATE_ASSOCIATION

#include "gin/stats/CChi2.h"
#include "gin/globals.h"

class UnivariateAssociation {
public:

	UnivariateAssociation();
	UnivariateAssociation(MatrixXd*, VectorXd*);

	// association scores
	virtual VectorXd computeSKAT();
	virtual VectorXd computeSKAT(VectorXd);
	virtual VectorXd computeChi2();
	virtual VectorXd computeTrendTest(std::string const&);

private:
	MatrixXd* __X;
	VectorXd* __y;
	uint __n_features;
	bool __binary_y;
};

#endif //GIN_UNIVARIATE_ASSOCIATION
