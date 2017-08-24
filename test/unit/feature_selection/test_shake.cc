//
// Created by hclimente on 31/07/2017.
//

#include "gin/feature_selection/shake.h"
#include "gin/gwas/CGWASData.h"
#include "gtest/gtest.h"

TEST(Shake, readGWAS) {

	Shake s;
	GWASData g;

	s.readGWAS("test/data/case1/genotype", 0);
	CPlinkParser::readPEDFile("test/data/case1/genotype.ped", &g);

	EXPECT_EQ(s.gwas()->raw_snps, g.raw_snps);
	EXPECT_EQ(s.gwas()->Y, g.Y);

	CPlinkParser::readMAPFile("test/data/case1/genotype.map", &g);

	EXPECT_EQ(s.gwas()->chromosomes, g.chromosomes);
	EXPECT_EQ(s.gwas()->positions, g.positions);
	EXPECT_EQ(s.gwas()->snp_identifiers, g.snp_identifiers);
	EXPECT_EQ(s.gwas()->snp_distance, g.snp_distance);

	CGWASDataHelper::encodeHeterozygousData(&g, 0);
	EXPECT_EQ(s.gwas()->X, g.X);

}

TEST(Shake, readNetwork) {

	Shake s;

	s.readGWAS("test/data/case1/genotype", 0);

	GWASData g = (*s.gwas());

	s.readNetwork("test/data/case1/network.txt");

	CSconesIO::readSparseNetworkFile("test/data/case1/network.txt", &g);
	EXPECT_EQ(MatrixXd(s.gwas()->network), MatrixXd(g.network));

}

// TODO selectHyperparameters, selectSNPs