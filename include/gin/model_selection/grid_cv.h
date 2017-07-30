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

	GridCV(MatrixXd* const&, VectorXd* const&, SparseMatrixXd* const&, VectorXd, uint, uint);
	GridCV(MatrixXd* const&, VectorXd* const&, SparseMatrixXd* const&, VectorXd, VectorXd, uint, uint);
	~GridCV();

	// grid exploration functions
	void exploreGrids(uint);
	double scoreModels(VectorXd const&, uint const&);

	// setters & getters
	std::pair<double, double> bestParameters() { return __bestParameters; }

private:

	// grid parameters
	std::vector<Grid*> __grids;
	uint __folds;
	VectorXd __etas;
	VectorXd __lambdas;

	// grid results
	std::map<double, std::map<double, VectorXd>> __gridAggregation;
	MatrixXd __gridSummary;
	std::pair<double, double> __bestParameters;

	// data
	MatrixXd* __X = NULL;
	VectorXd* __y = NULL;
	SparseMatrixXd* __W = NULL;
	bool __binary_y;

	// methods
	void __initGrids(uint);
	double __computeConsistency(VectorXd const&);
	double __computeInformation(VectorXd, uint);

};

#endif //GIN_GRID_CV_H_H
