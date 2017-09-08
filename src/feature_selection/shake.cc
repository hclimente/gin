//
// Created by hclimente on 21/07/2017.
//

#include "gin/feature_selection/shake.h"

Shake::Shake() {
	__gwas = new GWASData();
	__cvgrid = NULL;
	__debug = false;
}

Shake::Shake(MatrixXd X, VectorXd y, SparseMatrixXd W) {
	__gwas = new GWASData();
	__gwas->X = X;
	__gwas->y = y;
	__gwas->network = W;
	__cvgrid = NULL;
	__debug = false;
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

void Shake::selectHyperparameters(uint folds, uint const &modelScore, uint const &associationScore) {

	float64 begin = clock();
	__computeAssociation(associationScore);

	logging(STATUS,"Selecting the best hyperparameters...");
	__cvgrid = new GridCV(&(__gwas->X), &(__gwas->y), &(__gwas->network), folds);
	__cvgrid->initFolds(__c, associationScore);
	logging(INFO,"Running models.");
	__cvgrid->runFolds();
	logging(INFO,"Finding best model.");
	__cvgrid->scoreModels(modelScore);

	__bestEta = __cvgrid->bestEta();
	__bestLambda = __cvgrid->bestLambda();
	logging(WARNING,"Finished in " + StringHelper::to_string<float64>(float64(clock()-begin)/CLOCKS_PER_SEC) + " sec\n");

	if(__debug) {
		GridViews g(__cvgrid);
		logging(DEBUG, "Model score matrix");
		logging(DEBUG, g.viewSelectionCriterion());
		logging(DEBUG, "Average number of selected SNPs.");
		logging(DEBUG, g.viewSelectedAvg());
	}

}

void Shake::selectHyperparameters(uint folds, uint const &modelScore, uint const &associationScore, VectorXd const &etas, VectorXd const& lambdas) {

	float64 begin = clock();
	__computeAssociation(associationScore);

	logging(STATUS,"Selecting the best hyperparameters...");
	__cvgrid = new GridCV(&(__gwas->X), &(__gwas->y), &(__gwas->network), folds);
	__cvgrid->initFolds(etas, lambdas, associationScore);
	logging(INFO,"Running models.");
	__cvgrid->runFolds();
	logging(INFO,"Finding best model.");
	__cvgrid->scoreModels(modelScore);

	__bestEta = __cvgrid->bestEta();
	__bestLambda = __cvgrid->bestLambda();
	logging(WARNING,"Finished in " + StringHelper::to_string<float64>(float64(clock()-begin)/CLOCKS_PER_SEC) + " sec\n");

	if(__debug) {
		GridViews g(__cvgrid);
		logging(DEBUG, "Model score matrix");
		logging(DEBUG, g.viewSelectionCriterion());
		logging(DEBUG, "Average number of selected SNPs.");
		logging(DEBUG, g.viewSelectedAvg());
	}

}

void Shake::selectSNPs() {

	float64 begin = clock();
	logging(STATUS,"Searching ConES with eta = " + StringHelper::to_string<float64>(__bestEta) + " and lambda = " + StringHelper::to_string<float64>(__bestLambda) + "\n");
	Scones s = Scones(__c, __bestEta, __bestLambda, &(__gwas->network));
	s.selectSnps();
	__selectedSnps = s.selected();
	logging(WARNING,"Finished in " + StringHelper::to_string<float64>(float64(clock()-begin)/CLOCKS_PER_SEC) + " sec\n");

}

void Shake::writeResults(string const& output) {

	float64 begin = clock();
	logging(STATUS,"Writing results...");
	CSconesIO::writeOutput(output, __gwas, __selectedSnps, __bestEta, __bestLambda, __c);
	logging(WARNING,"Finished in " + StringHelper::to_string<float64>(float64(clock()-begin)/CLOCKS_PER_SEC) + " sec\n");

}

void Shake::__computeAssociation(uint const &associationScore) {

	logging(INFO,"Computing univariate association.");
	UnivariateAssociation univar( &(__gwas->X), &(__gwas->y) );

	if (associationScore == CHI2) {
		__c = univar.computeChi2();
	} else if (associationScore == TREND) {
		// TODO implement different trend scores
		__c = univar.computeTrendTest("additive");
	} else if (associationScore == SKAT) {
		__c = univar.computeSKAT();
	}

}