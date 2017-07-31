//
// Created by hclimente on 23/07/2017.
//

#ifndef GIN_GRID_CV_H_H
#define GIN_GRID_CV_H_H

#include "gin/model_selection/grid.h"
#include "gin/model_selection/CCrossValidation.h"
#include "gin/regression/CRegression.h"
#include "gin/utils/CMatrixHelper.h"

class GridCV {

public:

	GridCV(MatrixXd* const&, VectorXd* const&, SparseMatrixXd* const&, VectorXd, uint);
	GridCV(MatrixXd* const&, VectorXd* const&, SparseMatrixXd* const&, VectorXd, VectorXd, uint);
	~GridCV();

	// grid exploration functions
	void runFolds(uint);
	void scoreModels(uint);

	// setters & getters
	std::pair<double, double> bestParameters() { return __bestParameters; }
	VectorXd etas() { return __etas; }
	VectorXd lambdas() { return __lambdas; }
	std::vector<Grid*> grids() { return __grids; }
	bool binary_y() { return __binary_y; }
	std::map<double, std::map<double, VectorXd>> aggregatedFolds() { return __aggregatedFolds; };
	VectorXd aggregatedFolds(uint const& e, uint const& l) { return __aggregatedFolds[__etas(e)][__lambdas(l)]; };
	MatrixXd scoredFolds() { return __scoredFolds; }

	void set_binary_y(bool binary_y) { __binary_y = binary_y; }

private:

	// grid parameters
	std::vector<Grid*> __grids;
	uint __folds;
	VectorXd __etas;
	VectorXd __lambdas;

	// grid results
	std::map<double, std::map<double, VectorXd>> __aggregatedFolds;
	MatrixXd __scoredFolds;
	std::pair<double, double> __bestParameters;

	void __setAggregatedFolds(uint const &e, uint const &l, VectorXd const &val) { __aggregatedFolds[__etas(e)][__lambdas(l)] = val; }

	// data
	MatrixXd* __X;
	VectorXd* __y;
	SparseMatrixXd* __W;
	bool __binary_y;

	// methods
	double __computeConsistency(VectorXd const&);
	double __computeInformation(VectorXd, uint);

};

#endif //GIN_GRID_CV_H_H
