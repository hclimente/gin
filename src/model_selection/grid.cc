//
// Created by hclimente on 22/07/2017.
//

#include "gin/model_selection/grid.h"

Grid::Grid(MatrixXd const& X, VectorXd const& y, SparseMatrixXd* const& W ) {
	__X = X;
	__y = y;
	__W = W;

	__computeUnivariate(0);

	// TODO autoparams based on c
	__etas = VectorXd(10);
	__lambdas = VectorXd(10);

	__initGrid();

}

Grid::Grid(MatrixXd const& X, VectorXd const& y, VectorXd const& etas, VectorXd const& lambdas, SparseMatrixXd* const& W) {
	__X = X;
	__y = y;
	__etas = etas;
	__lambdas = lambdas;
	__W = W;
	__initGrid();
	__computeUnivariate(0);

}

void Grid::__initGrid() {

	for (int e = 0; e < __etas.rows(); e++) {
		for (int l = 0; l < __lambdas.rows(); l++) {
			double eta = __etas[e];
			double lambda = __lambdas[l];
			__grid[eta][lambda] = VectorXd();
		}
	}
}

void Grid::__computeUnivariate(uint const& association) {

	UnivariateAssociation ass(__X, __y);
	VectorXd c;

	if(association == SKAT) {
		__c = ass.computeSKAT();
	} else if (association == CHI2) {
		__c = ass.computeChi2();
	} else if (association == TREND) {
		__c = ass.computeTrendTest("additive");
	}

}

void Grid::search() {

	for (int e = 0; e < __etas.rows(); e++) {
		for (int l = 0; l < __lambdas.rows(); l++) {

			double eta = __etas[e];
			double lambda = __lambdas[l];

			Scones s(__c, eta, lambda, *__W);
			s.selectSnps();

			__grid[eta][lambda] = s.selected();

		}
	}
}

VectorXd Grid::selected(double const& eta, double const& lambda) {
	return __grid[eta][lambda];
}

std::vector<VectorXd> Grid::selected(VectorXd const &etas, VectorXd const &lambdas) {
	std::vector<VectorXd> selected;

	for (int e = 0; e < etas.rows(); e++) {
		for (int l = 0; l < lambdas.rows(); l++) {
			double eta = etas[e];
			double lambda = lambdas[l];
			selected.push_back(__grid[eta][lambda]);
		}
	}

	return selected;
}