//
// Created by hclimente on 22/07/2017.
//

#include "gin/model_selection/grid.h"

Grid::Grid(VectorXd const& c, SparseMatrixXd* const& W) {

	__c = c;
	__W = W;


	__etas = VectorXd::LinSpaced(10, log10(__c.maxCoeff()), log10(__c.minCoeff()));
	for(int i = 0; i < __etas.rows(); i++)
		__etas(i) = pow(10, __etas(i));

	__lambdas = VectorXd::LinSpaced(10, log10(__c.maxCoeff()), log10(__c.minCoeff()));
	for(int i = 0; i < __lambdas.rows(); i++)
		__lambdas(i) = pow(10, __lambdas(i));

	__initGrid();

}

Grid::Grid(VectorXd const& c, SparseMatrixXd* const& W, VectorXd const& etas, VectorXd const& lambdas) {

	__c = c;
	__etas = etas;
	__lambdas = lambdas;
	__W = W;
	__initGrid();

	// sort etas decreasingly
	std::sort(__etas.data(), __etas.data() + __etas.size(), [](int a, int b){return a > b;});

}

void Grid::__initGrid() {

	for (uint64 l = 0; l < __lambdas.rows(); l++) {
		for (uint64 e = 0; e < __etas.rows(); e++) {
			double eta = __etas[e];
			double lambda = __lambdas[l];
			__grid[eta][lambda] = VectorXd();
		}
	}
}

bool Grid::__trivial(double const& lambda) {

	VectorXd trivial = VectorXd();

	// look for trivial cases in the solution boundaries
	// case with most sparsity (largest eta): lower bound for the solution
	double upper_eta = __etas(0);
	Scones upper_bound(__c, upper_eta, lambda, __W);
	upper_bound.selectSnps();
	__grid[upper_eta][lambda] = upper_bound.selected();

	// case with least sparsity (smallest eta): upper bound for the solution
	double lower_eta = __etas(__etas.rows() - 1);
	Scones lower_bound(__c, lower_eta, lambda, __W);
	lower_bound.selectSnps();
	__grid[lower_eta][lambda] = lower_bound.selected();

	// check for trivial solutions before calculating
	if (__grid[lower_eta][lambda] == VectorXd::Zero(__c.rows())) {
		trivial = __grid[lower_eta][lambda];
	} else if (__grid[upper_eta][lambda] == VectorXd::Ones(__c.rows())) {
		trivial = __grid[upper_eta][lambda];
	}

	for (uint64 e = 1; e < __etas.rows() - 1; e++) {
		double eta = __etas(e);
		__grid[eta][lambda] = trivial;
	}

	return (trivial != VectorXd());

}

void Grid::search() {

	for (uint64 l = 0; l < __lambdas.rows(); l++) {

		double lambda = __lambdas(l);

		if (__trivial(lambda)) {
			continue;
		}

		for (uint64 e = 1; e < __etas.rows() - 1; e++) {

			double eta = __etas[e];

			Scones s(__c, eta, lambda, __W);
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

	for (uint64 l = 0; l < lambdas.rows(); l++) {
		for (uint64 e = 0; e < etas.rows(); e++) {
			double eta = etas[e];
			double lambda = lambdas[l];
			selected.push_back(__grid[eta][lambda]);
		}
	}

	return selected;
}