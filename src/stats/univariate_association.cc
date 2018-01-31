//
// Created by hclimente on 21/07/2017.
//

#include "gin/stats/univariate_association.h"
#include "gin/regression/CRegression.h"

UnivariateAssociation::UnivariateAssociation() {
	__X = NULL;
	__y = NULL;
	__n_features = 0;

	__binary_y = true;

}

UnivariateAssociation::UnivariateAssociation(MatrixXd* X, VectorXd* y) {
	__X = X;
	__y = y;
	__n_features = __X -> cols();

	__binary_y = true;
	for(int64 i = 0; i < __y -> rows(); i++) {
		// check if 0=unknown, 1=unaffected or 2=affected
		if(!( (*__y)(i) == 0 || (*__y)(i) == 1 || (*__y)(i) == 2)) {
			__binary_y = false;
			break;
		}
	}
}

VectorXd UnivariateAssociation::computeSKAT() {

	VectorXd r = __y->array() - __y->mean();

	VectorXd nonweighted_skat = ((__X->transpose()*r).array().pow(2));
	return nonweighted_skat;

}

VectorXd UnivariateAssociation::computeSKAT(VectorXd W) {

	MatrixXd sW;
	sW = DiagXd(__n_features);
	sW.diagonal() = W;

	VectorXd nonweighted_skat = computeSKAT();
	return sW * nonweighted_skat;
}

VectorXd UnivariateAssociation::computeChi2() {

	if (!__binary_y) {
		// TODO throw error
	}

	VectorXd chi2(__n_features);

	for (uint64 i = 0; i < __n_features; i++) {
		MatrixXd tab = CChi2::get2DContingencyTable(__X -> col(i), (*__y), true );
		chi2(i) = CChi2::calculateChi2(tab);
	}

	return chi2;

}

VectorXd UnivariateAssociation::computeTrendTest(std::string const& geneticModel) {

	if (!__binary_y) {
		// TODO throw error
	}

	VectorXd T(__n_features);
	VectorXd model(3);

	if(geneticModel == "dominance")
		model << 1, 1, 0;
	else if(geneticModel == "recessive")
		model << 0, 1, 1;
	else if(geneticModel == "codominance")
		model << 0, 1, 2;
	else
		model << 1, 1, 1;

	for (uint64 i = 0; i < __X -> cols(); i++) {
		MatrixXd tab = CChi2::get2DContingencyTable(__X -> col(i), (*__y), true);
		T(i) = CChi2::calculateChi2Trend(tab, model);
	}

	return T;
}