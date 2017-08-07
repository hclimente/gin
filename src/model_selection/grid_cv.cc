//
// Created by hclimente on 24/07/2017.
//

#include "gin/model_selection/grid_cv.h"

GridCV::GridCV(MatrixXd* const& X, VectorXd* const& y, SparseMatrixXd* const& W, VectorXd c, uint folds) {
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
		// check if 0=unknown, 1=unaffected or 2=affected
		if(!( (*__y)(i) == 0 || (*__y)(i) == 1 || (*__y)(i) == 2)) {
			__binary_y = false;
			break;
		}
	}

	if(__binary_y) {
		// TODO check how to make logistic regression less computationally expensive
		// __regressor = new CLogisticRegression();
		__regressor = new CLinearRegression();
	} else {
		__regressor = new CLinearRegression();
	}

	__folds = folds;
	__scoredFolds = MatrixXd::Zero(__etas.rows(), __lambdas.rows());

}

GridCV::GridCV(MatrixXd* const& X, VectorXd* const& y, SparseMatrixXd* const& W, VectorXd etas, VectorXd lambdas, uint folds) {
	__X = X;
	__y = y;
	__W = W;

	__binary_y = true;
	for(int64 i = 0; i < __y->rows(); i++) {
		// check if 0=unknown, 1=unaffected or 2=affected
		if(!( (*__y)(i) == 0 || (*__y)(i) == 1 || (*__y)(i) == 2)) {
			__binary_y = false;
			break;
		}
	}

	if(__binary_y) {
		// TODO check how to make logistic regression less computationally expensive
		// __regressor = new CLogisticRegression();
		__regressor = new CLinearRegression();
	} else {
		__regressor = new CLinearRegression();
	}

	__etas = etas;
	__lambdas = lambdas;
	__folds = folds;
	__scoredFolds = MatrixXd::Zero(__etas.rows(), __lambdas.rows());

}

GridCV::~GridCV() {

	for (int i = 0; i < __grids.size(); i++) {
		delete __grids[i];
	}

	delete __regressor;

}

void GridCV::runFolds(uint association) {

	// delete previous iteration, if any
	for (int i = 0; i < __grids.size(); i++) {
		delete __grids[i];
		__grids[i] = NULL;
	}

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
			__setAggregatedFolds(e, l, VectorXd());
			__scoredFolds(e,l) = -1;
		}
	}

	for (int i = 0; i < __folds; i++) {
		__grids[i] -> search();
	}
}

void GridCV::scoreModels(uint scoring_function) {

	for (uint e = 0; e < __etas.rows(); e++ ) {
		for (uint l = 0; l < __lambdas.rows(); l++ ) {

			__setAggregatedFolds(e, l, VectorXd::Zero(__X->cols()));

			for (int i = 0; i < __folds; i++) {
				VectorXd u = __grids[i] -> selected(__etas[e], __lambdas[l]);
				__setAggregatedFolds(e, l, aggregatedFolds(e, l) + u);
			}
			__setAggregatedFolds(e, l, aggregatedFolds(e, l) / __folds);

			// TODO consider cases where a number > that a threshold are picked
			if (scoring_function == CONSISTENCY) {
				__scoredFolds(e, l) = __computeConsistency(aggregatedFolds(e,l));
			} else if(scoring_function == AICc | scoring_function == AIC | scoring_function == BIC | scoring_function == mBIC) {
				__scoredFolds(e, l) = - __computeInformation(aggregatedFolds(e,l), scoring_function);
			}
		}
	}

	MatrixXd::Index best_eta_index, best_lambda_index;
	__scoredFolds.maxCoeff(&best_eta_index, &best_lambda_index);
	__bestEta = __etas[best_eta_index];
	__bestLambda = __lambdas[best_lambda_index];

}

double GridCV::__computeConsistency(VectorXd const& folds) {

	// TODO implementation is not the same as in original SConES
	double score;
	float totalSelected = (folds.array() > 0).count() * __folds;

	if (totalSelected == 0 ) {
		score = 0;
	} else {
		float consistentlySelected = (folds.array() * __folds).sum();
		score = consistentlySelected / totalSelected;
	}

	return score;

}

double GridCV::__computeInformation(VectorXd folds, uint scoringFunction) {

    // max possible value if no features selected
	double score = 1e31;

	for(uint64 i = 0; i < folds.rows(); i++) {
		if(folds(i) < 1) {
			folds(i) = 0;
		}
	}
	
	if (folds.sum() != 0) {
		MatrixXd x_tr = sliceColsMatrixByBinaryVector(*__X, folds);

		__regressor->fit(*__y, x_tr);
		
		if (scoringFunction == BIC) {
			score = __regressor->getBIC();
		} else if(scoringFunction == AIC) {
			score = __regressor->getAIC();
		}  else if(scoringFunction == AICc) {
			score = __regressor->getAICc();
		}  else if(scoringFunction == mBIC) {
			// TODO implement network information
			score = __regressor->getBIC();
		}
	}

	return score;

}
