//
// Created by hclimente on 23/08/2017.
//

#ifndef GIN_GRID_VIEWS_H
#define GIN_GRID_VIEWS_H

#include "gin/model_selection/grid_cv.h"
#include "gin/globals.h"

class GridViews {

friend GridCV;

public:

	GridViews(GridCV* const&);
	~GridViews();

	MatrixXd viewSelectionCriterion();
	MatrixXd viewSelectedAvg();

private:

	GridCV* __grid_cv;

	MatrixXd __baseView();

};

#endif //GIN_GRID_VIEWS_H
