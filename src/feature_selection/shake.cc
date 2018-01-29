//
// Created by hclimente on 21/07/2017.
//

#include "gin/feature_selection/shake.h"

Shake::Shake() {
	__gwas = new GWASData();
	__cvgrid = NULL;
	__debug = false;
}

Shake::Shake(MatrixXd* X, VectorXd* y, SparseMatrixXd* W) {
	__gwas = new GWASData();
	__gwas->X = *X;
	__gwas->y = *y;
	__gwas->network = *W;
	__cvgrid = NULL;
	__debug = false;
}

Shake::~Shake() {
	delete __gwas;
	delete __cvgrid;
}

void Shake::readGWAS(string const& pedBasename, uint encoding) {

	float64 begin = clock();
	logging(GIN_STATUS,"Reading Genotype file...");
	CPlinkParser::readPEDFile(pedBasename + ".ped", __gwas);
	logging(GIN_INFO,"Number of SNPs: " + StringHelper::to_string<uint64>(__gwas->n_snps));
	logging(GIN_INFO,"Number of Samples: " + StringHelper::to_string<uint64>(__gwas->n_samples));
	logging(GIN_WARNING,"Finished in " + StringHelper::to_string<float64>(float64(clock()-begin)/CLOCKS_PER_SEC) + " sec\n");

	begin = clock();
	logging(GIN_STATUS,"Reading Mapping file...");
	CPlinkParser::readMAPFile(pedBasename + ".map", __gwas);
	logging(GIN_WARNING,"Finished in " + StringHelper::to_string<float64>(float64(clock()-begin)/CLOCKS_PER_SEC) + " sec\n");

	begin = clock();
	logging(GIN_STATUS,"Encoding SNP data...");
	CGWASDataHelper::encodeHeterozygousData(__gwas, encoding);
	logging(GIN_WARNING,"Finished in " + StringHelper::to_string<float64>(float64(clock()-begin)/CLOCKS_PER_SEC) + " sec\n");

}

void Shake::readNetwork(string const& networkFilename) {

	float64 begin = clock();
	logging(GIN_STATUS,"Loading and filtering network file...");
	CSconesIO::readSparseNetworkFile(networkFilename, __gwas);
	logging(GIN_WARNING,"Finished in " + StringHelper::to_string<float64>(float64(clock()-begin)/CLOCKS_PER_SEC) + " sec\n");

}

void Shake::selectHyperparameters(uint folds, uint const &modelScore, uint const &associationScore) {

	float64 begin = clock();
	__computeAssociation(associationScore);

	logging(GIN_STATUS,"Selecting the best hyperparameters...");
	__cvgrid = new GridCV(&(__gwas->X), &(__gwas->y), &(__gwas->network), folds);
	__cvgrid->initFolds(__c, associationScore);
	logging(GIN_INFO,"Running models.");
	__cvgrid->runFolds();
	logging(GIN_INFO,"Finding best model.");
	__cvgrid->scoreModels(modelScore);

	__bestEta = __cvgrid->bestEta();
	__bestLambda = __cvgrid->bestLambda();
	logging(GIN_WARNING,"Finished in " + StringHelper::to_string<float64>(float64(clock()-begin)/CLOCKS_PER_SEC) + " sec\n");

	if(__debug) {
		GridViews g(__cvgrid);
		logging(GIN_DEBUG, "Model score matrix");
		logging(GIN_DEBUG, g.viewSelectionCriterion());
		logging(GIN_DEBUG, "Average number of selected SNPs.");
		logging(GIN_DEBUG, g.viewSelectedAvg());
	}

}

void Shake::selectHyperparameters(uint folds, uint const &modelScore, uint const &associationScore, VectorXd const &etas, VectorXd const& lambdas) {

	float64 begin = clock();
	__computeAssociation(associationScore);

	logging(GIN_STATUS,"Selecting the best hyperparameters...");
	__cvgrid = new GridCV(&(__gwas->X), &(__gwas->y), &(__gwas->network), folds);
	__cvgrid->initFolds(etas, lambdas, associationScore);
	logging(GIN_INFO,"Running models.");
	__cvgrid->runFolds();
	logging(GIN_INFO,"Finding best model.");
	__cvgrid->scoreModels(modelScore);

	__bestEta = __cvgrid->bestEta();
	__bestLambda = __cvgrid->bestLambda();
	logging(GIN_WARNING,"Finished in " + StringHelper::to_string<float64>(float64(clock()-begin)/CLOCKS_PER_SEC) + " sec\n");

	if(__debug) {
		GridViews g(__cvgrid);
		logging(GIN_DEBUG, "Model score matrix");
		logging(GIN_DEBUG, g.viewSelectionCriterion());
		logging(GIN_DEBUG, "Average number of selected SNPs.");
		logging(GIN_DEBUG, g.viewSelectedAvg());
	}

}

void Shake::selectSNPs() {

	float64 begin = clock();
	logging(GIN_STATUS,"Searching ConES with eta = " + StringHelper::to_string<float64>(__bestEta) + " and lambda = " + StringHelper::to_string<float64>(__bestLambda) + "\n");
	Scones s = Scones(__c, __bestEta, __bestLambda, &(__gwas->network));
	s.selectSnps();
	__selectedSnps = s.selected();
	logging(GIN_WARNING,"Finished in " + StringHelper::to_string<float64>(float64(clock()-begin)/CLOCKS_PER_SEC) + " sec\n");

}

void Shake::writeResults(string const& output) {

	float64 begin = clock();
	logging(GIN_STATUS,"Writing results...");
	CSconesIO::writeOutput(output, __gwas, __selectedSnps, __bestEta, __bestLambda, __c);
	logging(GIN_WARNING,"Finished in " + StringHelper::to_string<float64>(float64(clock()-begin)/CLOCKS_PER_SEC) + " sec\n");

}

void Shake::__computeAssociation(uint const &associationScore) {

	logging(GIN_INFO,"Computing univariate association.");
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