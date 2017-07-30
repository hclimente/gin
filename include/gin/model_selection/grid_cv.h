//
// Created by hclimente on 23/07/2017.
//

#ifndef GIN_GRID_CV_H_H
#define GIN_GRID_CV_H_H

#include "gin/model_selection/grid.h"
#include "gin/model_selection/CCrossValidation.h"
#include "gin/utils/CMatrixHelper.h"

class GridCV {

public:

	GridCV(MatrixXd* const&, VectorXd* const&, SparseMatrixXd* const&, VectorXd, uint);
	GridCV(MatrixXd* const&, VectorXd* const&, SparseMatrixXd* const&, VectorXd, VectorXd, uint);
	~GridCV();

	// grid exploration functions
	void exploreGrids(std::string const&);
	double scoreModels(VectorXd const&, std::string const&);

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

	// methods
	void __initGrids();

};

#endif //GIN_GRID_CV_H_H
