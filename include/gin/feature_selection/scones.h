//
// Created by hclimente on 21/07/2017.
//

#ifndef GIN_SCONES_H
#define GIN_SCONES_H

#include "gin/globals.h"
#include "gin/feature_selection/feature_selector.h"

class Scones: public FeatureSelector {
public:

	Scones(VectorXd const&, double const&, double const&, SparseMatrixXd* const&);
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
