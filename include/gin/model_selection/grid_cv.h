//
// Created by hclimente on 23/07/2017.
//

#ifndef GIN_GRID_CV_H
#define GIN_GRID_CV_H

#include "gin/model_selection/grid.h"
#include "gin/model_selection/CCrossValidation.h"
#include "gin/regression/CRegression.h"
#include "gin/stats/univariate_association.h"
#include "gin/utils/CMatrixHelper.h"

class GridCV {

public:

	GridCV();
	GridCV(MatrixXd* const&, VectorXd* const&, SparseMatrixXd* const&, uint);
	virtual ~GridCV();

	// grid exploration functions
	virtual void initFolds(VectorXd, VectorXd, uint);
	virtual void initFolds(VectorXd, uint);
	virtual void initFolds(uint);
	virtual void runFolds();
	virtual void scoreModels(uint);

	// setters & getters
	virtual double bestEta() { return __bestEta; }
	virtual double bestLambda() { return __bestLambda; }
	virtual VectorXd etas() { return __etas; }
	virtual VectorXd lambdas() { return __lambdas; }
	virtual std::vector<Grid*> grids() { return __grids; }
	virtual bool binary_y() { return __binary_y; }
	virtual std::map<double, std::map<double, VectorXd> > aggregatedFolds() { return __aggregatedFolds; };
	virtual VectorXd aggregatedFolds(uint const& e, uint const& l) { return __aggregatedFolds[__etas(e)][__lambdas(l)]; };
	virtual MatrixXd scoredFolds() { return __scoredFolds; }
	virtual double scoredFolds(uint const& e, uint const& l) { return __scoredFolds(e,l); }

	virtual void setCV(CCrossValidation cv) { __cv = cv; }

	void set_binary_y(bool binary_y) { __binary_y = binary_y; }

private:

	// grid parameters
	std::vector<Grid*> __grids;
	uint __folds;
	VectorXd __etas;
	VectorXd __lambdas;

	// grid results
	std::map<double, std::map<double, VectorXd> > __aggregatedFolds;
	MatrixXd __scoredFolds;
	double __bestEta;
	double __bestLambda;

	void __setAggregatedFolds(uint const &e, uint const &l, VectorXd const &val) { __aggregatedFolds[__etas(e)][__lambdas(l)] = val; }

	// data
	MatrixXd* __X;
	VectorXd* __y;
	SparseMatrixXd* __W;
	bool __binary_y;
	CRegression* __classifier;
	CCrossValidation __cv;

	// methods
	double __computeConsistency(VectorXd const&);
	double __computeInformation(VectorXd, uint);
	VectorXd __computeUnivariate(MatrixXd* const&, VectorXd* const&, uint const&);
	void __getClassifier();

};

#endif //GIN_GRID_CV_H_H
