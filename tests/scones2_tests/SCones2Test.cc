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

};

struct SearchMarkers : public ::testing::Test, testing::WithParamInterface<CSconesInitialSettings> {

    CSconesSettings settings;
    CScones* scones;
    float64 eta;
    float64 lambda;

    SearchMarkers(){

        GWASData data;

        string genotype_str = "data/testing/scones/genotype";
        string phenotype_str = "data/testing/scones/phenotype.txt";
        string network_str = "data/testing/scones/network.txt";
        uint encoding = 0;
        float64 maf = 0.05;

        CPlinkParser::readPEDFile(genotype_str + ".ped", &data);
        CPlinkParser::readMAPFile(genotype_str + ".map", &data);
        CPlinkParser::readPhenotypeFile(phenotype_str,&data);
        CGWASDataHelper::encodeHeterozygousData(&data,encoding);
        CGWASDataHelper::filterSNPsByMAF(&data,maf);
        CSconesIO::readSparseNetworkFile(network_str,&data);
        GWASData tmpData = CGWASDataHelper::removeSamples4MissingData(data,0);

        settings = CSconesSettings();
        eta = GetParam().eta;
        lambda = GetParam().lambda;

        settings.folds = 10;
        settings.seed = 0;
        settings.selection_criterion = CONSISTENCY;
        settings.selection_ratio = 0.8;
        settings.test_statistic = SKAT;
        settings.nParameters = 10;
        settings.evaluateObjective = false;
        settings.dump_intermediate_results = true;
        settings.dump_path = "tmp/";

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

struct CSconesTest : public ::testing::Test {
    CSconesSettings settings;
    CScones* scones;

    CSconesTest()
    {
        GWASData data;

        string genotype_str = "data/testing/scones/genotype";
        string phenotype_str = "data/testing/scones/phenotype.txt";
        string network_str = "data/testing/scones/network.txt";
        uint encoding = 0;
        float64 maf = 0.05;

        CPlinkParser::readPEDFile(genotype_str + ".ped", &data);
        CPlinkParser::readMAPFile(genotype_str + ".map", &data);
        CPlinkParser::readPhenotypeFile(phenotype_str,&data);
        CGWASDataHelper::encodeHeterozygousData(&data,encoding);
        CGWASDataHelper::filterSNPsByMAF(&data,maf);
        CSconesIO::readSparseNetworkFile(network_str,&data);
        GWASData tmpData = CGWASDataHelper::removeSamples4MissingData(data,0);

        settings = CSconesSettings();
        settings.folds = 10;
        settings.seed = 0;
        settings.selection_criterion = CONSISTENCY;
        settings.selection_ratio = 0.8;
        settings.test_statistic = SKAT;
        settings.autoParameters = true;
        settings.nParameters = 10;
        settings.lambdas = VectorXd::Zero(settings.nParameters);
        settings.etas = VectorXd::Zero(settings.nParameters);
        settings.evaluateObjective = false;
        settings.dump_intermediate_results = true;
        settings.dump_path = "tmp/";

        scones = new CScones(tmpData.Y.col(0),tmpData.X,tmpData.network, settings);
    }

    ~CSconesTest()
    {
        delete scones;
    }
};

struct CSconesTest_FixedLambdaEta : public ::testing::Test {
    CSconesSettings settings;
    CScones* scones;

    CSconesTest_FixedLambdaEta()
    {
        GWASData data;

        string genotype_str = "data/testing/scones/genotype";
        string phenotype_str = "data/testing/scones/phenotype.txt";
        string network_str = "data/testing/scones/network.txt";
        uint encoding = 0;
        float64 maf = 0.05;

        CPlinkParser::readPEDFile(genotype_str + ".ped", &data);
        CPlinkParser::readMAPFile(genotype_str + ".map", &data);
        CPlinkParser::readPhenotypeFile(phenotype_str,&data);
        CGWASDataHelper::encodeHeterozygousData(&data,encoding);
        CGWASDataHelper::filterSNPsByMAF(&data,maf);
        CSconesIO::readSparseNetworkFile(network_str,&data);
        GWASData tmpData = CGWASDataHelper::removeSamples4MissingData(data,0);

        settings = CSconesSettings();
        settings.folds = 10;
        settings.seed = 0;
        settings.selection_criterion = CONSISTENCY;
        settings.selection_ratio = 0.8;
        settings.test_statistic = SKAT;
        settings.nParameters = 10;
        settings.evaluateObjective = false;
        settings.dump_intermediate_results = true;
        settings.dump_path = "tmp/";

        // set up specific lambda and eta
        VectorXd l(1);
        l(0) = 300;
        VectorXd e(1);
        e(0) = 17000;
        settings.lambdas = l;
        settings.etas = e;
        // avoid gridsearch
        settings.autoParameters = false;

        scones = new CScones(tmpData.Y.col(0),tmpData.X,tmpData.network, settings);
    }

    ~CSconesTest_FixedLambdaEta()
    {
        delete scones;
    }
};

TEST_F(CSconesTest, integrity_selectedSNPs) {
    scones -> test_associations();
    int out = scones -> getIndicatorVector().sum();
    EXPECT_EQ(10, out);
    EXPECT_EQ(1, scones -> getIndicatorVector()(676));
    EXPECT_EQ(1, scones -> getIndicatorVector()(679));
    EXPECT_EQ(1, scones -> getIndicatorVector()(680));
    EXPECT_EQ(1, scones -> getIndicatorVector()(682));
    EXPECT_EQ(1, scones -> getIndicatorVector()(684));
    EXPECT_EQ(1, scones -> getIndicatorVector()(685));
    EXPECT_EQ(1, scones -> getIndicatorVector()(686));
    EXPECT_EQ(1, scones -> getIndicatorVector()(690));
    EXPECT_EQ(1, scones -> getIndicatorVector()(695));
    EXPECT_EQ(1, scones -> getIndicatorVector()(696));

}

TEST_F(CSconesTest, integrity_bestLambdaAndEta) {
    scones -> test_associations();
    float64 lambda = scones -> getBestLambda();
    float64 eta = scones -> getBestEta();

    EXPECT_NEAR(278.25, lambda , 0.1);
    EXPECT_NEAR(16681, eta, 0.1);

}

TEST_F(CSconesTest, integrity_objectiveFunctionTerms) {
    scones -> test_associations();
    float64 lambda = scones -> getBestLambda();
    float64 eta = scones -> getBestEta();
    VectorXd terms = scones -> getObjectiveFunctionTerms(lambda,eta);

    double association = terms(0);
    double connectivity = terms(1);
    double sparsity = terms(2);

    EXPECT_NEAR(724379, association, 1);
    EXPECT_NEAR(44799.20, connectivity, 1);
    EXPECT_NEAR(166810, sparsity, 1);

}

TEST_F(CSconesTest_FixedLambdaEta, integrity_bestLambdaAndEta_fixedLambdaAndEta) {
    scones -> test_associations();
    float64 lambda = scones -> getBestLambda();
    float64 eta = scones -> getBestEta();

    EXPECT_NEAR(300, lambda , 1);
    EXPECT_NEAR(17000, eta, 1);

}

TEST_P(SearchMarkers, check_lambda_eta) {
    auto as = GetParam();

    scones -> test_associations();
    float64 lambda = scones -> getBestLambda();
    float64 eta = scones -> getBestEta();

    EXPECT_NEAR(as.expected_lambda, lambda , 1);
    EXPECT_NEAR(as.expected_eta, eta, 1);
}

TEST_P(SearchMarkers, checkMarkers) {
    scones -> test_associations();
    int out = scones -> getIndicatorVector().sum();
    EXPECT_EQ(10, out);
    EXPECT_EQ(1, scones -> getIndicatorVector()(676));
    EXPECT_EQ(1, scones -> getIndicatorVector()(679));
    EXPECT_EQ(1, scones -> getIndicatorVector()(680));
    EXPECT_EQ(1, scones -> getIndicatorVector()(682));
    EXPECT_EQ(1, scones -> getIndicatorVector()(684));
    EXPECT_EQ(1, scones -> getIndicatorVector()(685));
    EXPECT_EQ(1, scones -> getIndicatorVector()(686));
    EXPECT_EQ(1, scones -> getIndicatorVector()(690));
    EXPECT_EQ(1, scones -> getIndicatorVector()(695));
    EXPECT_EQ(1, scones -> getIndicatorVector()(696));
}

INSTANTIATE_TEST_CASE_P(checkParameters, SearchMarkers,
    testing::Values(
            CSconesInitialSettings{-1, 278.25, -1, 16681},
            CSconesInitialSettings{300, 300, 17000, 17000}
    ));


int main(int argc, char* argv[])
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}