//
// Created by hclimente on 23/08/2017.
//

#include "gin/io/grid_views.h"

GridViews::GridViews(GridCV* const& gridcv) {

	__grid_cv = gridcv;

}

GridViews::~GridViews() {

}

MatrixXd GridViews::viewSelectionCriterion() {

	MatrixXd view = __baseView();

	for (uint e = 0; e < __grid_cv->etas().rows(); e++ ) {
		for (uint l = 0; l < __grid_cv->lambdas().rows(); l++ ) {
			view(e+1,l+1) = __grid_cv->scoredFolds(e,l);
		}
	}

	return view;

}

MatrixXd GridViews::viewSelectedAvg() {

	MatrixXd view = __baseView();

	for (uint e = 0; e < __grid_cv->etas().rows(); e++ ) {
		for (uint l = 0; l < __grid_cv->lambdas().rows(); l++ ) {
			view(e+1,l+1) = __grid_cv->aggregatedFolds(e,l).sum();
		}
	}

	return view;

}

MatrixXd GridViews::__baseView() {

	MatrixXd baseView = MatrixXd::Zero(__grid_cv->etas().rows() + 1, __grid_cv->lambdas().rows() + 1);
	baseView.block(1, 0, __grid_cv->etas().rows(), 1) = __grid_cv->etas();
	baseView.block(0, 1, 1, __grid_cv->lambdas().rows()) = __grid_cv->lambdas().transpose();

	return baseView;

}