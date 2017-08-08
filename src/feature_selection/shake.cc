//
// Created by hclimente on 21/07/2017.
//

#include "gin/feature_selection/shake.h"

Shake::Shake() {
	__gwas = new GWASData();
	__cvgrid = NULL;
}

Shake::Shake(MatrixXd X, VectorXd Y, SparseMatrixXd W) {
	__gwas = new GWASData();
	__gwas->X = X;
	__gwas->Y = Y;
	__gwas->network = W;
	__cvgrid = NULL;
}

Shake::~Shake() {
	delete __gwas;
	delete __cvgrid;
}

void Shake::readGWAS(string const& pedBasename, uint encoding) {

	float64 begin = clock();
	logging(STATUS,"Reading Genotype file...");
	CPlinkParser::readPEDFile(pedBasename + ".ped", __gwas);
	logging(INFO,"Number of SNPs: " + StringHelper::to_string<uint64>(__gwas->n_snps));
	logging(INFO,"Number of Samples: " + StringHelper::to_string<uint64>(__gwas->n_samples));
	logging(WARNING,"Finished in " + StringHelper::to_string<float64>(float64(clock()-begin)/CLOCKS_PER_SEC) + " sec\n");

	begin = clock();
	logging(STATUS,"Reading Mapping file...");
	CPlinkParser::readMAPFile(pedBasename + ".map", __gwas);
	logging(WARNING,"Finished in " + StringHelper::to_string<float64>(float64(clock()-begin)/CLOCKS_PER_SEC) + " sec\n");

	begin = clock();
	logging(STATUS,"Encoding SNP data...");
	CGWASDataHelper::encodeHeterozygousData(__gwas, encoding);
	logging(WARNING,"Finished in " + StringHelper::to_string<float64>(float64(clock()-begin)/CLOCKS_PER_SEC) + " sec\n");

}

void Shake::readNetwork(string const& networkFilename) {

	float64 begin = clock();
	logging(STATUS,"Loading and filtering network file...");
	CSconesIO::readSparseNetworkFile(networkFilename, __gwas);
	logging(WARNING,"Finished in " + StringHelper::to_string<float64>(float64(clock()-begin)/CLOCKS_PER_SEC) + " sec\n");

}

void Shake::searchHyperparameters(uint folds, uint const& modelScore, uint const& associationScore) {

	float64 begin = clock();
	logging(STATUS,"Selecting the best hyperparameters...");
	MatrixXd X = __gwas -> X;
	VectorXd y = __gwas -> Y.col(0);
	SparseMatrixXd W = __gwas -> network;

	logging(INFO,"Computing univariate association.");
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
	logging(INFO,"Running models.");
	__cvgrid->runFolds(associationScore);
	logging(INFO,"Finding best models.");
	__cvgrid->scoreModels(modelScore);

	__bestEta = __cvgrid->bestEta();
	__bestLambda = __cvgrid->bestLambda();
	logging(WARNING,"Finished in " + StringHelper::to_string<float64>(float64(clock()-begin)/CLOCKS_PER_SEC) + " sec\n");

}

void Shake::selectSnps() {

	float64 begin = clock();
	logging(STATUS,"Searching ConES...\n");
	SparseMatrixXd W = __gwas -> network;

	Scones s = Scones(__c, __bestEta, __bestLambda, &W);
	__selectedSnps = s.selected();
	logging(WARNING,"Finished in " + StringHelper::to_string<float64>(float64(clock()-begin)/CLOCKS_PER_SEC) + " sec\n");

}

void Shake::writeResults(string const& output) {

	float64 begin = clock();
	logging(STATUS,"Writing results...");
	CSconesIO::writeOutput(output, __gwas, __selectedSnps, __bestEta, __bestLambda, __c);
	logging(WARNING,"Finished in " + StringHelper::to_string<float64>(float64(clock()-begin)/CLOCKS_PER_SEC) + " sec\n");

}