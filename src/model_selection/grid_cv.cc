//
// Created by hclimente on 24/07/2017.
//

#include "gin/model_selection/grid_cv.h"

GridCV::GridCV() {
	__X = NULL;
	__y = NULL;
	__W = NULL;

	__folds = 0;
	__etas = VectorXd();
	__lambdas = VectorXd();

	__classifier = NULL;
	__scoredFolds = MatrixXd::Zero(__etas.rows(), __lambdas.rows());
}

GridCV::GridCV(MatrixXd* const& X, VectorXd* const& y, SparseMatrixXd* const& W, uint folds) {

	// constructor just sets up the data, and partitions it in the k-fold setting
	__X = X;
	__y = y;
	__W = W;
	__folds = folds;

	// TODO get seed from user settings
	__cv = CCrossValidation(0);
	__cv.kFold(__folds, __y->rows());

	__getClassifier();

}

GridCV::~GridCV() {

	for (int i = 0; i < __grids.size(); i++) {
		delete __grids[i];
	}

	delete __classifier;

}

void GridCV::initFolds(VectorXd etas, VectorXd lambdas, uint associationScore) {

	__etas = etas;
	__lambdas = lambdas;
	__scoredFolds = MatrixXd::Zero(__etas.rows(), __lambdas.rows());

	initFolds(associationScore);

}

void GridCV::initFolds(VectorXd c, uint associationScore) {

	double minc = c.minCoeff();

	// change by small quantity to avoid problems with log(0)
	if (minc == 0) {
		minc = minc + 0.000001;
	}

	double maxc = c.maxCoeff();

	__etas = VectorXd::LinSpaced(10, log10(minc), log10(maxc));
	for(int i = 0; i < __etas.rows(); i++)
		__etas(i) = pow(10, __etas(i));

	__lambdas = VectorXd::LinSpaced(10, log10(minc) - 1, log10(maxc) + 1);
	for(int i = 0; i < __lambdas.rows(); i++)
		__lambdas(i) = pow(10, __lambdas(i));

	__scoredFolds = MatrixXd::Zero(__etas.rows(), __lambdas.rows());

	initFolds(associationScore);

}

void GridCV::initFolds(uint associationScore) {

	// delete previous iteration, if any
	for (int i = 0; i < __grids.size(); i++) {
		delete __grids[i];
		__grids[i] = NULL;
	}

	for (int i = 0; i < __folds; i++) {
		VectorXd tr_indices = __cv.getTrainingIndices(i);
		MatrixXd X_train = sliceRowsMatrix(*__X, tr_indices);
		VectorXd y_train = sliceRowsMatrix(*__y, tr_indices);

		VectorXd c = __computeUnivariate(&X_train, &y_train, associationScore);
		/* TODO what is this
		MatrixXd cov_train;
		if(__covs_set) {
			cov_train = sliceRowsMatrix(__covs,tr_indices);
		}
		 */

		Grid* g = new Grid(c, __W, __etas, __lambdas);
		__grids.push_back(g);
	}

}

void GridCV::runFolds() {

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

VectorXd GridCV::__computeUnivariate(MatrixXd* const& X, VectorXd* const& y, uint const& association) {

	UnivariateAssociation snpAssociation(X, y);
	VectorXd c;

	if(association == SKAT) {
		c = snpAssociation.computeSKAT();
	} else if (association == CHI2) {
		c = snpAssociation.computeChi2();
	} else if (association == TREND) {
		c = snpAssociation.computeTrendTest("additive");
	}

	return c;

}