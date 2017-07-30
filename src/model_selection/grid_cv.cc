//
// Created by hclimente on 24/07/2017.
//

#include "gin/model_selection/grid_cv.h"

GridCV::GridCV(MatrixXd* const& X, VectorXd* const& y, SparseMatrixXd* const& W, VectorXd c, uint folds) {
	__X = X;
	__y = y;
	__W = W;

	// TODO extend range
	// autoparams
	__etas = VectorXd::LinSpaced(10, log10(c.maxCoeff()), log10(c.minCoeff()));
	for(int i = 0; i < __etas.rows(); i++)
		__etas(i) = pow(10, __etas(i));

	__lambdas = VectorXd::LinSpaced(10, log10(c.maxCoeff()), log10(c.minCoeff()));
	for(int i = 0; i < __lambdas.rows(); i++)
		__lambdas(i) = pow(10, __lambdas(i));

	__folds = folds;
	__gridSummary = MatrixXd::Zero(__etas.rows(), __lambdas.rows());
	__initGrids();
}

GridCV::GridCV(MatrixXd* const& X, VectorXd* const& y, SparseMatrixXd* const& W, VectorXd etas, VectorXd lambdas, uint folds) {
	__X = X;
	__y = y;
	__W = W;
	__etas = etas;
	__lambdas = lambdas;
	__folds = folds;
	__gridSummary = MatrixXd::Zero(__etas.rows(), __lambdas.rows());
	__initGrids();
}

GridCV::~GridCV() {

	for (int i = 0; i < __folds; i++) {
		delete __grids[i];
	}

}

void GridCV::__initGrids() {

	// TODO get seed
	CCrossValidation cv(0);
	cv.kFold(__folds, __y -> rows());

	for (int i = 0; i < __folds; i++) {
		VectorXd tr_indices = cv.getTrainingIndices(i);
		MatrixXd x_train = sliceRowsMatrix(*__X, tr_indices);
		VectorXd y_train = sliceRowsMatrix(*__y, tr_indices);
		/* TODO what is this
		MatrixXd cov_train;
		if(__covs_set) {
			cov_train = sliceRowsMatrix(__covs,tr_indices);
		}
		 */

		Grid* g = new Grid(x_train, y_train, __etas, __lambdas, __W);
		__grids.push_back(g);
	}

	for (int e = 0; e < __etas.rows(); e++ ) {
		for (int l = 0; l < __lambdas.rows(); l++ ) {
			__gridAggregation[e][l] = VectorXd();
			__gridSummary(e, l) = -1;
		}
	}
}

void GridCV::exploreGrids(std::string const& scoring_function) {

	for (int i = 0; i < __folds; i++) {
		__grids[i] -> search();
	}

	double max = -1;

	for (int e = 0; e < __etas.rows(); e++ ) {
		for (int l = 0; l < __lambdas.rows(); l++ ) {

			double eta = __etas[e];
			double lambda = __lambdas[e];
			__gridAggregation[eta][lambda] = VectorXd::Zero(__X -> cols());

			for (int i = 0; i < __folds; i++) {
				VectorXd u = __grids[i] -> getSelected(eta, lambda);
				__gridAggregation[eta][lambda] += u;
			}

			// TODO implement the other scoring functions
			__gridSummary(e, l) = scoreModels(__gridAggregation[eta][lambda], scoring_function);

			if(max < __gridSummary(e, l)) {
				max = __gridSummary(e, l);
				__bestParameters = std::pair<double, double>(eta, lambda);
			}
		}
	}
}

double GridCV::scoreModels(VectorXd const& folds, std::string const& scoringFunction) {

	double score;

	// TODO implementation is not the same as in original SConES
	if (scoringFunction == "CONSISTENCY") {
		float totalSelected = (folds.array() > 0).count();

		if (totalSelected == 0 ) {
			score = 0;
		} else {
			float consistentlySelected = (folds.array() == __folds).count();
			score = consistentlySelected / totalSelected;
		}
	}

	return score;

}