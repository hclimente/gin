#include "gtest/gtest.h"
#include "CEasyGWAS/gwas/CScones.h"
#include "CEasyGWAS/io/CSconesIO.h"
#include "CEasyGWAS/io/CPlinkParser.h"
#include "CEasyGWAS/globals.h"

struct CSconesInitialSettings{

    float64 eta;
    float64 expected_eta;
    float64 lambda;
    float64 expected_lambda;
    int griddepth;
    uint test_statistic;
    float64 expected_association;
    float64 expected_connectivity;
    float64 expected_sparsity;
    string path_prefix;
    int selection_criterion;
    vector<int> expected_causal_SNPs;

};

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

CSconesInitialSettings skat_grid_consistency = CSconesInitialSettings {
        -1, // eta
        16681.00537, // expected_eta
        -1, // lambda
        278.25594, // expected_lambda
        1, // griddepth
        SKAT, // test_statistic
        724379, // expected_association
        44799.20, // expected_connectivity
        166810, // expected_sparsity
        "data/testing/scones/skat/", // path_prefix
        CONSISTENCY, // selection_criterion
        {676, 679, 680, 682, 684, 685, 686, 690, 695, 696} // expected_causal_SNPs
};

CSconesInitialSettings skat_grid_information = CSconesInitialSettings {
        -1, // eta
        16681.00537, // expected_eta
        -1, // lambda
        2154.43469, // expected_lambda
        1, // griddepth
        SKAT, // test_statistic
        820863, // expected_association
        4308, // expected_connectivity
        433706, // expected_sparsity
        "data/testing/scones/skat/", // path_prefix
        AICc, // selection_criterion
        {676, 677, 678, 679, 680, 681, 682, 683, 684, 685, 686, 687, 688, 689, 690, 691, 692, 693, 694, 695, 696, 697, 698, 699, 700, 701} // expected_causal_SNPs
};

CSconesInitialSettings skat_fixed = CSconesInitialSettings {
        17000, // eta
        17000, // expected_eta
        300, // lambda
        300, // expected_lambda
        1, // griddepth
        SKAT, // test_statistic
        724379, // expected_association
        48300, // expected_connectivity
        170000, // expected_sparsity
        "data/testing/scones/skat/", // path_prefix
        CONSISTENCY, // selection_criterion
        {676, 679, 680, 682, 684, 685, 686, 690, 695, 696} // expected_causal_SNPs
};

CSconesInitialSettings chisq_fixed = CSconesInitialSettings {
        13.0694, // eta
        13.0694, // expected_eta
        0.218, // lambda
        0.218, // expected_lambda
        1, // griddepth
        CHISQ, // test_statistic
        517.83, // expected_association
        35.1, // expected_connectivity
        130.694, // expected_sparsity
        "data/testing/scones/skat/", // path_prefix
        CONSISTENCY, // selection_criterion
        {676, 679, 680, 682, 684, 685, 686, 690, 695, 696} // expected_causal_SNPs

};

CSconesInitialSettings chisq_grid5_information = CSconesInitialSettings {
        -1, // eta
        33435.468, // expected_eta
        -1, // lambda
        3.385, // expected_lambda
        5, // griddepth
        SKAT, // test_statistic
        625692, // expected_association
        453.695, // expected_connectivity
        234048, // expected_sparsity
        "data/testing/scones/skat/", // path_prefix
        CONSISTENCY, // selection_criterion
        {676, 680, 682, 685, 686, 690, 696} // expected_causal_SNPs

};

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