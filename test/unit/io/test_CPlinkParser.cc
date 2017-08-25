//
// Created by hclimente on 01/08/2017.
//

#include "gin/io/CPlinkParser.h"
#include "gin/gwas/CGWASData.h"
#include "gtest/gtest.h"
#include "gin/globals.h"


TEST(CPlinkParser, readPEDFile) {

	GWASData g;
	CPlinkParser::readPEDFile("test/data/test_chunks/read.ped", &g);

	EXPECT_EQ(g.n_snps, 3);
	EXPECT_EQ(g.n_samples, 10);
	EXPECT_EQ(g.raw_snps[0][0], 'C');
	EXPECT_EQ(g.raw_snps[0][1], 'M');
	EXPECT_EQ(g.raw_snps[1][0], 'M');
	EXPECT_EQ(g.raw_snps[1][1], 'C');
	EXPECT_EQ(g.raw_snps[2][0], 'Y');
	EXPECT_EQ(g.raw_snps[2][1], 'Y');
	EXPECT_EQ(g.raw_snps[3][0], 'S');
	EXPECT_EQ(g.raw_snps[3][1], 'S');
	EXPECT_EQ(g.raw_snps[4][0], 'A');
	EXPECT_EQ(g.raw_snps[4][1], 'W');
	EXPECT_EQ(g.raw_snps[5][0], 'W');
	EXPECT_EQ(g.raw_snps[5][1], 'A');
	EXPECT_EQ(g.raw_snps[6][0], 'R');
	EXPECT_EQ(g.raw_snps[6][1], 'T');
	EXPECT_EQ(g.raw_snps[7][0], 'T');
	EXPECT_EQ(g.raw_snps[7][1], 'R');
	EXPECT_EQ(g.raw_snps[8][0], 'G');
	EXPECT_EQ(g.raw_snps[8][1], 'K');
	EXPECT_EQ(g.raw_snps[9][0], 'K');
	EXPECT_EQ(g.raw_snps[9][1], 'G');

	VectorXd y(10);
	y << .1, .2, .3, .4, .5, .6, .7, .8, .9, 1;

	for (uint i = 0; i < g.n_samples; i++) {
		EXPECT_NEAR(g.y(i), y[i], 0.01);
	}
}

