#include "gtest/gtest.h"
#include "CEasyGWAS/gwas/CScones.h"
#include "CEasyGWAS/io/CSconesIO.h"
#include "CEasyGWAS/io/CPlinkParser.h"
#include "CEasyGWAS/globals.h"

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

TEST_F(CSconesTest, integrity_objectiveFunctionTerms) {
    scones -> test_associations();
    float64 lambda = scones -> getBestLambda();
    float64 eta = scones -> getBestEta();
    VectorXd terms = scones -> getObjectiveFunctionTerms(lambda,eta);

    double association = terms(0);
    double connectivity = terms(1);
    double sparsity = terms(2);

    EXPECT_NEAR(724379, association, 1);
    EXPECT_NEAR(25043, connectivity, 1);
    EXPECT_NEAR(166810, sparsity, 1);

}

int main(int argc, char* argv[])
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}