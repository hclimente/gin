//
// Created by hclimente on 22/01/2018.
//

#include "gin/globals.h"

#ifndef AS_RGINLIB

void set_seed(float64 seed) {
	srand(seed);
}

VectorXd shuffle_vector(int size) {

	VectorXd indices = VectorXd::LinSpaced(size, 0, size - 1);
	Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> perm(size);
	perm.setIdentity();
	std::random_shuffle(perm.indices().data(), perm.indices().data()+perm.indices().size());
	indices = perm * indices;

	return indices;

}

int urand() {

	int x = rand();
	return x;

}

VectorXd urand(int size) {

	VectorXd x = VectorXd::Random(size);
	return x;

}

#else
#include "Rcpp.h"
#include <math.h>

void set_seed(float64 seed) {}

inline int randWrapper(const int n) { return floor(unif_rand()*n); }

VectorXd shuffle_vector(int size) {

	VectorXd indices = VectorXd::LinSpaced(size, 0, size - 1);
	Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> perm(size);
	perm.setIdentity();
	std::random_shuffle(perm.indices().data(), perm.indices().data()+perm.indices().size(), randWrapper);
	indices = perm * indices;

	return indices;

}

int urand() {

	int x = floor(R::runif(0, RAND_MAX));
	return x;

}

VectorXd urand(int size) {

	VectorXd x = VectorXd::Zero(size);
	for (int i = 0; i < size; i++) {
		x(i) = R::runif(-1,1);
	}
	return x;

}
#endif //AS_RGINLIB
