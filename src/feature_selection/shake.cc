//
// Created by hclimente on 21/07/2017.
//

#include "gin/feature_selection/shake.h"

Shake::Shake() {
	__gwas = new GWASData();
	__cvgrid = NULL;
}

Shake::Shake(GWASData* gwas) {
	__gwas = gwas;
	__cvgrid = NULL;
}

Shake::~Shake() {
	delete __gwas;
	delete __cvgrid;
}

void Shake::readGWAS(string const& pedBasename, uint encoding) {

	CPlinkParser::readPEDFile(pedBasename + ".ped", __gwas);
	CPlinkParser::readMAPFile(pedBasename + ".map", __gwas);

	CGWASDataHelper::encodeHeterozygousData(__gwas, encoding);
}

void Shake::readNetwork(string const& networkFilename) {

	CSconesIO::readSparseNetworkFile(networkFilename, __gwas);

}

void Shake::searchHyperparameters(uint folds, uint const& modelScore, uint const& associationScore) {

	MatrixXd X = __gwas -> X;
	VectorXd y = __gwas -> Y.col(0);
	SparseMatrixXd W = __gwas -> network;

	UnivariateAssociation univar( &X, &y );

	if (associationScore == CHI2) {
		__c = univar.computeChi2();
	} else if (associationScore == TREND) {
		// TODO implement different trend scores
		__c = univar.computeTrendTest("additive");
	} else if (associationScore == SKAT) {
		__c = univar.computeSKAT();
	}

	__cvgrid = new GridCV(&X, &y, &W, __c, folds);
	__cvgrid->runFolds(associationScore);
	__cvgrid->scoreModels(modelScore);

	__bestEta = __cvgrid -> bestParameters().first;
	__bestLambda = __cvgrid -> bestParameters().second;

}

void Shake::selectSnps() {

	SparseMatrixXd W = __gwas -> network;

	Scones s = Scones(__c, __bestEta, __bestLambda, &W);
	__selectedSnps = s.selected();

}

void Shake::writeResults(string const& output) {

	CSconesIO::writeOutput(output, __gwas, __selectedSnps, __bestEta, __bestLambda, __c);

}