//
// Created by hclimente on 21/07/2017.
//

#ifndef GIN_SCONES_H
#define GIN_SCONES_H

#include "gin/globals.h"
#include "gin/feature_selection/feature_selector.h"

class Scones: public FeatureSelector {
public:
	Scones(VectorXd const& c, double const& eta, double const& lambda, SparseMatrixXd const& W)
			: FeatureSelector::FeatureSelector(c.rows())
	{

		// TODO throw exception in f and nrow do not match
		__c = c;
		__eta = eta;
		__lambda = lambda;
		__lW = __lambda * W;
	};
	void selectSnps();
	double computeScore();

private:
	VectorXd __c;
	double __eta;
	double __lambda;
	SparseMatrixXd __lW;

	void __maxflow(MatrixXd const &);

};

#endif //GIN_SCONES_H
