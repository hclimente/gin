//
// Created by hclimente on 21/07/2017.
//

#include "gin/feature_selection/scones.h"
#include "maxflow/maxflow.h"

Scones::Scones(VectorXd const& c, double const& eta, double const& lambda, SparseMatrixXd* const& W)
: FeatureSelector::FeatureSelector(c.rows())
{

	if (c.rows() != W -> rows()){
		// TODO throw exception in f and nrow do not match
	}

	__c = c;
	__eta = eta;
	__lambda = lambda;
	__lW = __lambda * (*W);

}

void Scones::selectSnps() {

	VectorXd c_t = __c.array() - __eta;

	// Add source and sink
	// Matrix A containing As and At
	MatrixXd A(__n_features, 2);
	// connect positive c values to sink
	VectorXd pos_c = (c_t.array() <= 0).select(0, c_t);
	// connect negative c values to source
	VectorXd neg_c = - c_t;
	neg_c = (c_t.array() > 0).select(0, neg_c);
	//Store data
	A.col(0) = neg_c;
	A.col(1) = pos_c;

	//compute maxflow
	__maxflow(A);

}

void Scones::__maxflow(MatrixXd const &A) {

	// create graph out of adjacency matrix
	typedef Graph<float64, float64, float64> MaxGraph;
	MaxGraph *g = new MaxGraph(__lW.rows(), __lW.nonZeros());

	// initialize nodes
	g->add_node(__selected.rows());

	//traverse the sparse adjacency matrix A
	for(int64 k=0; k<__lW.outerSize(); k++) {
		for(Eigen::SparseMatrix<float64>::InnerIterator it(__lW,k); it; ++it) {
			g->add_edge(it.row(), it.col(), it.value(), 0.0);
		}
	}

	// traverse the T matrix
	for(int64 i = 0; i < 2; i++) {
		for(int64 k = 0; k < A.rows(); k++) {
			if(i==0) {
				if(A(k,i) != 0) {
					g->add_tweights(k, A(k,i), 0.0);
				}
			} else {
				if(A(k,i) != 0) {
					g->add_tweights(k, 0.0, A(k,i));
				}
			}
		}
	}

	// run maxflow algorithm
	g->maxflow();

	// create indicator_vector
	for(int64 i = 0; i < __selected.rows(); i++)
		__selected(i) = g->what_segment(i);

	// delete graph
	delete g;

}

double Scones::computeScore() {

	double score = 0.0;

	for(uint i = 0; i < __c.rows(); i++) {
		if(__selected(i) == 1) {
			score += __c(i);
			score -= __eta;
		}
	}

	for(int64 k = 0; k < __lW.outerSize(); k++) {
		for(Eigen::SparseMatrix<float64>::InnerIterator it(__lW, k); it; ++it) {
			if(__selected(it.row()) == 1) {
				if(__selected(it.col()) == 0) {
					score -= it.value();
				}
			}
		}
	}

	return score;
}