//
// Created by hclimente on 22/01/2018.
//

#include "gin/globals.h"

#ifndef AS_RGINLIB

void set_seed(float64 seed) {
	srand(seed);
}

VectorXd shuffle_vector(int N) {

	VectorXd indices = VectorXd::LinSpaced(N, 0, N - 1);
	Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> perm(N);
	perm.setIdentity();
	std::random_shuffle(perm.indices().data(), perm.indices().data()+perm.indices().size());
	indices = perm * indices;

	return indices;

}

#else
#include "Rcpp.h"

void set_seed(float64 seed) {}

inline int randWrapper(const int n) { return floor(unif_rand()*n); }

VectorXd shuffle_vector(int N) {

	VectorXd indices = VectorXd::LinSpaced(N, 0, N - 1);
	Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> perm(N);
	perm.setIdentity();
	std::random_shuffle(perm.indices().data(), perm.indices().data()+perm.indices().size(), randWrapper);
	indices = perm * indices;

	return indices;

}
#endif //AS_RGINLIB

VectorXd shuffle_vector(VectorXd);
