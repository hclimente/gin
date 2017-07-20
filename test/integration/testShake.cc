#include "testShake.h"

struct SearchMarkers : public ::testing::Test, testing::WithParamInterface<CSconesInitialSettings> {

    CSconesSettings settings;
    CScones* scones;
    float64 eta;
    float64 lambda;
    GWASData tmpData;

    SearchMarkers(){

        GWASData data;

        string genotype_str = GetParam().path_prefix + "genotype";
        string phenotype_str = GetParam().path_prefix + "phenotype.txt";
        string network_str = GetParam().path_prefix + "network.txt";
        uint encoding = 0;
        float64 maf = 0.05;

        CPlinkParser::readPEDFile(genotype_str + ".ped", &data);
        CPlinkParser::readMAPFile(genotype_str + ".map", &data);
        CPlinkParser::readPhenotypeFile(phenotype_str,&data);
        CGWASDataHelper::encodeHeterozygousData(&data,encoding);
        CGWASDataHelper::filterSNPsByMAF(&data,maf);
        CSconesIO::readSparseNetworkFile(network_str,&data);
        tmpData = CGWASDataHelper::removeSamples4MissingData(data,0);

        eta = GetParam().eta;
        lambda = GetParam().lambda;

        settings = CSconesSettings();
        settings.folds = 10;
        settings.seed = 0;
        settings.selection_criterion = GetParam().selection_criterion;
        settings.selection_ratio = 0.8;
        settings.test_statistic = GetParam().test_statistic;
        settings.nParameters = 10;
        settings.evaluateObjective = true;
        settings.dump_intermediate_results = true;
        settings.dump_path = "tmp/";
        settings.gridsearch_depth = GetParam().griddepth;

        if (eta >= 0 & lambda >= 0) {
            // set up specific lambda and eta
            VectorXd l(1);
            l(0) = lambda;
            VectorXd e(1);
            e(0) = eta;
            settings.lambdas = l;
            settings.etas = e;
            // avoid gridsearch
            settings.autoParameters = false;
        } else {
            settings.lambdas = VectorXd::Zero(settings.nParameters);
            settings.etas = VectorXd::Zero(settings.nParameters);
            settings.autoParameters = true;
        }

        scones = new CScones(tmpData.Y.col(0),tmpData.X,tmpData.network, settings);
    }

    ~SearchMarkers(){
        delete scones;
    }
};

TEST_P(SearchMarkers, checkObjectiveFunctionTerms) {
    auto as = GetParam();
    scones -> test_associations();
    float64 lambda = scones -> getBestLambda();
    float64 eta = scones -> getBestEta();
    VectorXd terms = scones -> getObjectiveFunctionTerms(lambda,eta);

    double association = terms(0);
    double connectivity = terms(1);
    double sparsity = terms(2);

    EXPECT_NEAR(as.expected_association, association, 1);
    EXPECT_NEAR(as.expected_connectivity, connectivity, 1);
    EXPECT_NEAR(as.expected_sparsity, sparsity, 1);

}

TEST_P(SearchMarkers, checkLambdaAndEta) {
    auto as = GetParam();

    scones -> test_associations();
    float64 lambda = scones -> getBestLambda();
    float64 eta = scones -> getBestEta();

    EXPECT_NEAR(as.expected_lambda, lambda , 0.001);
    EXPECT_NEAR(as.expected_eta, eta, 0.001);
}

TEST_P(SearchMarkers, checkSelectedSNPS) {
    auto as = GetParam();

    scones -> test_associations();
    int out = scones -> getIndicatorVector().sum();
    EXPECT_EQ(as.expected_causal_SNPs.size(), out);
    for (unsigned int i = 0; i < as.expected_causal_SNPs.size(); i++)
        EXPECT_EQ(1, scones -> getIndicatorVector()(as.expected_causal_SNPs[i]));
}

TEST_P(SearchMarkers, checkSelectedSNPs_fixedParameters) {
    auto as = GetParam();

    scones -> test_associations(as.expected_lambda, as.expected_eta);
    VectorXd indicator = scones -> getIndicatorVector();

    EXPECT_EQ(as.expected_causal_SNPs.size(), indicator.sum());
    for (unsigned int i = 0; i < as.expected_causal_SNPs.size(); i++){
        EXPECT_EQ(1, indicator(as.expected_causal_SNPs[i]));
    }
}

TEST_P(SearchMarkers, checkOutputFiles) {
    auto as = GetParam();

    scones -> test_associations();
    VectorXd terms = scones -> getObjectiveFunctionTerms(scones -> getBestLambda(), scones -> getBestEta());
    VectorXd skat = scones -> getScoreStatistic();

    string output_str = as.path_prefix + "/" + tmpData.phenotype_names[0] + ".scones.out.ext.txt";
    CSconesIO::writeOutput(output_str, tmpData, scones -> getIndicatorVector(), scones -> getBestLambda(), scones -> getBestEta(), terms, skat);

    output_str = as.path_prefix  + "/" + tmpData.phenotype_names[0] + ".scones.pmatrix.txt";
    CSconesIO::writeCMatrix(output_str, scones -> getCMatrix(), scones -> getSettings());
}

INSTANTIATE_TEST_CASE_P(checkParameters, SearchMarkers,
    testing::Values(
            skat_fixed,
            skat_grid_consistency,
            skat_grid_information,
            chisq_fixed,
            chisq_grid5_information
    ));


int main(int argc, char* argv[])
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}