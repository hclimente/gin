//
// Created by hclimente on 24/07/2017.
//

#include "gin/model_selection/grid_cv.h"

GridCV::GridCV() {
	__X = NULL;
	__y = NULL;
	__W = NULL;

	__folds = 0;
	__etas = VectorXd::Zero(1);
	__lambdas = VectorXd::Zero(1);

	__classifier = NULL;
	__scoredFolds = MatrixXd::Zero(__etas.rows(), __lambdas.rows());
}

GridCV::GridCV(MatrixXd* const& X, VectorXd* const& y, SparseMatrixXd* const& W, VectorXd c, uint folds) {

	__X = X;
	__y = y;
	__W = W;

	__etas = VectorXd::LinSpaced(10, log10(c.minCoeff()), log10(c.maxCoeff()));
	for(int i = 0; i < __etas.rows(); i++)
		__etas(i) = pow(10, __etas(i));

	__lambdas = VectorXd::LinSpaced(10, log10(c.minCoeff()) - 1, log10(c.maxCoeff()) + 1);
	for(int i = 0; i < __lambdas.rows(); i++)
		__lambdas(i) = pow(10, __lambdas(i));

	__getClassifier();

	__folds = folds;
	__scoredFolds = MatrixXd::Zero(__etas.rows(), __lambdas.rows());

}

GridCV::GridCV(MatrixXd* const& X, VectorXd* const& y, SparseMatrixXd* const& W, VectorXd etas, VectorXd lambdas, uint folds) {

	__X = X;
	__y = y;
	__W = W;

	__getClassifier();

	__etas = etas;
	__lambdas = lambdas;
	__folds = folds;
	__scoredFolds = MatrixXd::Zero(__etas.rows(), __lambdas.rows());

}

GridCV::~GridCV() {

	for (int i = 0; i < __grids.size(); i++) {
		delete __grids[i];
	}

	delete __classifier;

}

void GridCV::runFolds(uint association) {

	// TODO get seed from user settings
	CCrossValidation cv(0);
	cv.kFold(__folds, __y->rows());
	runFolds(association, cv);

}

void GridCV::runFolds(uint association, CCrossValidation cv) {

	// delete previous iteration, if any
	for (int i = 0; i < __grids.size(); i++) {
		delete __grids[i];
		__grids[i] = NULL;
	}

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
				__setAggregatedFolds(e, l, aggregatedFolds(e,l) + u);
			}
			__setAggregatedFolds(e, l, aggregatedFolds(e,l) / __folds);

			if (aggregatedFolds(e,l).array().sum() > (__X->cols() * 0.05) | aggregatedFolds(e,l).array().sum() == 0) {
				__scoredFolds(e, l) = -1e31;
			} else {
				if (scoring_function == CONSISTENCY) {
					__scoredFolds(e, l) = __computeConsistency(aggregatedFolds(e, l));
				} else if (scoring_function == AICc | scoring_function == AIC | scoring_function == BIC |
				           scoring_function == mBIC) {
					__scoredFolds(e, l) = -__computeInformation(aggregatedFolds(e, l), scoring_function);
				}
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

		__classifier->fit(*__y, x_tr);
		
		if (scoringFunction == BIC) {
			score = __classifier->getBIC();
		} else if(scoringFunction == AIC) {
			score = __classifier->getAIC();
		}  else if(scoringFunction == AICc) {
			score = __classifier->getAICc();
		}  else if(scoringFunction == mBIC) {
			// TODO implement network information
			score = __classifier->getBIC();
		}
	}

	return score;

}

void GridCV::__getClassifier() {

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
		// __classifier = new CLogisticRegression();
		__classifier = new CLinearRegression();
	} else {
		__classifier = new CLinearRegression();
	}
}