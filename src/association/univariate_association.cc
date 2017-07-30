//
// Created by hclimente on 21/07/2017.
//

#include "gin/association/univariate_association.h"
#include "gin/regression/CRegression.h"

UnivariateAssociation::UnivariateAssociation(MatrixXd X, VectorXd y) {
	__X = X;
	__y = y;
	__n_features = __X.cols();

	__binary_y = true;
	for(int64 i = 0; i < __y.rows(); i++) {
		if(!(__y(i) == 0 || __y(i) == 1)) {
			__binary_y = false;
			break;
		}
	}
}

VectorXd UnivariateAssociation::computeSKAT() {

	MatrixXd sW;
	sW = DiagXd(__n_features);
	sW.diagonal() = VectorXd::Ones(__n_features);

	VectorXd nonweighted_skat = ((__X.transpose()*__y).array().pow(2));
	return sW * nonweighted_skat;
}

VectorXd UnivariateAssociation::computeChi2() {

	VectorXd chi2(__n_features);

	for (int i = 0; i < __n_features; i++) {
		MatrixXd tab = CChi2::get2DContingencyTable(__X.col(i), __y);
		chi2(i) = CChi2::calculateChi2(tab);
	}

	return chi2;

}

VectorXd UnivariateAssociation::computeTrendTest(std::string const& geneticModel) {
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

	for (int i = 0; i < __X.cols(); i++) {
		MatrixXd tab = CChi2::get2DContingencyTable(__y, __X.col(i));
		T(i) = CChi2::calculateChi2Trend(tab, model);
	}

	return T;
}