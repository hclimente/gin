//
// Created by hclimente on 24/07/2017.
//

#include "gin/model_selection/grid_cv.h"

GridCV::GridCV(MatrixXd* const& X, VectorXd* const& y, SparseMatrixXd* const& W, VectorXd c, uint folds, uint association) {
	__X = X;
	__y = y;
	__W = W;

	__etas = VectorXd::LinSpaced(10, log10(c.maxCoeff()), log10(c.minCoeff()));
	for(int i = 0; i < __etas.rows(); i++)
		__etas(i) = pow(10, __etas(i));

	// TODO extend range for lambdas
	__lambdas = VectorXd::LinSpaced(10, log10(c.maxCoeff()), log10(c.minCoeff()));
	for(int i = 0; i < __lambdas.rows(); i++)
		__lambdas(i) = pow(10, __lambdas(i));

	__binary_y = true;
	for(int64 i = 0; i < __y->rows(); i++) {
		if( !( (*__y)(i) == 0 || (*__y)(i) == 1)) {
			__binary_y = false;
			break;
		}
	}

	__folds = folds;
	__gridSummary = MatrixXd::Zero(__etas.rows(), __lambdas.rows());
	__initGrids(association);
}

GridCV::GridCV(MatrixXd* const& X, VectorXd* const& y, SparseMatrixXd* const& W, VectorXd etas, VectorXd lambdas, uint folds, uint association) {
	__X = X;
	__y = y;
	__W = W;

	__binary_y = true;
	for(int64 i = 0; i < __y->rows(); i++) {
		if( !( (*__y)(i) == 0 || (*__y)(i) == 1)) {
			__binary_y = false;
			break;
		}
	}

	__etas = etas;
	__lambdas = lambdas;
	__folds = folds;
	__gridSummary = MatrixXd::Zero(__etas.rows(), __lambdas.rows());
	__initGrids(association);
}

GridCV::~GridCV() {

	for (int i = 0; i < __folds; i++) {
		delete __grids[i];
	}

}

void GridCV::__initGrids(uint association) {

	// TODO get seed from user settings
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

		Grid* g = new Grid(x_train, y_train, __W, association, __etas, __lambdas);
		__grids.push_back(g);
	}

	for (int e = 0; e < __etas.rows(); e++ ) {
		for (int l = 0; l < __lambdas.rows(); l++ ) {
			__gridAggregation[e][l] = VectorXd();
			__gridSummary(e, l) = -1;
		}
	}
}

void GridCV::exploreGrids(uint scoring_function) {

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
				VectorXd u = __grids[i]->selected(eta, lambda);
				__gridAggregation[eta][lambda] += u;
			}

			__gridAggregation[eta][lambda] /= __folds;
			__gridSummary(e, l) = scoreModels(__gridAggregation[eta][lambda], scoring_function);

			if(max < __gridSummary(e, l)) {
				max = __gridSummary(e, l);
				__bestParameters = std::pair<double, double>(eta, lambda);
			}
		}
	}
}

double GridCV::scoreModels(VectorXd const& folds, uint const& scoringFunction) {

	double score;

	// TODO consider cases where a number > that a threshold are picked
	if (scoringFunction == CONSISTENCY) {
		score = __computeConsistency(folds);
	} else if(scoringFunction == AICc | scoringFunction == AIC | scoringFunction == BIC | scoringFunction == mBIC) {
		score = __computeInformation(folds, scoringFunction);
	}

	return score;

}

double GridCV::__computeConsistency(VectorXd const& folds) {

	// TODO implementation is not the same as in original SConES
	double score;
	float totalSelected = (folds.array() > 0).count();

	if (totalSelected == 0 ) {
		score = 0;
	} else {
		float consistentlySelected = (folds.array() == 1).count();
		score = consistentlySelected / totalSelected;
	}

	return score;

}

double GridCV::__computeInformation(VectorXd folds, uint scoringFunction) {

	double score;

	for(uint64 i = 0; i < folds.rows(); i++) {
		if(folds(i) != 1) {
			folds(i) = 0;
		}
	}

	// max possible value if no features selected
	if (folds.sum() == 0)
		score = 1e31;
	else {
		MatrixXd x_tr = sliceColsMatrixByBinaryVector(*__X, folds);
		CRegression regression;

		if(__binary_y) {
			regression = CLogisticRegression();
		} else {
			regression = CLinearRegression();
		}

		regression.fit(*__y, x_tr);

		if (scoringFunction == BIC) {
			score = regression.getBIC();
		} else if(scoringFunction == AIC) {
			score = regression.getAIC();
		}  else if(scoringFunction == AICc) {
			score = regression.getAICc();
		}  else if(scoringFunction == mBIC) {
			// TODO implement network information
			score = regression.getBIC();
		}
	}

	return score;

}