//
// Created by hclimente on 02/08/2017.
//

#include "gin/gwas/CGWASData.h"
#include "gin/io/CPlinkParser.h"
#include "gtest/gtest.h"
#include "gin/globals.h"

TEST(CGWASDataHelper, encodeHeterozygousData) {

	GWASData g;
	CPlinkParser::readPEDFile("test/data/test_chunks/encoding.ped", &g);

	// additive: each allele has a contribution to the phenotype
	CGWASDataHelper::encodeHeterozygousData(&g, 0);

	EXPECT_EQ(g.X(0,0), 0);
	EXPECT_EQ(g.X(1,0), 0);
	EXPECT_EQ(g.X(2,0), 1);
	EXPECT_EQ(g.X(3,0), 2);
	EXPECT_EQ(g.X(0,1), 2);
	EXPECT_EQ(g.X(1,1), 1);
	EXPECT_EQ(g.X(2,1), 0);
	EXPECT_EQ(g.X(3,1), 0);
	// no basis por picking one or the other in the last SNP
	EXPECT_EQ(g.X(0,2), 2);
	EXPECT_EQ(g.X(1,2), 2);
	EXPECT_EQ(g.X(2,2), 0);
	EXPECT_EQ(g.X(3,2), 0);

	// recessive: the minor allele causes a phenotype when in homozygosis
	CGWASDataHelper::encodeHeterozygousData(&g, 1);

	EXPECT_EQ(g.X(0,0), 0);
	EXPECT_EQ(g.X(1,0), 0);
	EXPECT_EQ(g.X(2,0), 0);
	EXPECT_EQ(g.X(3,0), 1);
	EXPECT_EQ(g.X(0,1), 1);
	EXPECT_EQ(g.X(1,1), 0);
	EXPECT_EQ(g.X(2,1), 0);
	EXPECT_EQ(g.X(3,1), 0);
	// no basis por picking one or the other in the last SNP
	EXPECT_EQ(g.X(0,2), 1);
	EXPECT_EQ(g.X(1,2), 1);
	EXPECT_EQ(g.X(2,2), 0);
	EXPECT_EQ(g.X(3,2), 0);

	// dominant: the minor allele causes a phenotype when present
	CGWASDataHelper::encodeHeterozygousData(&g, 2);

	EXPECT_EQ(g.X(0,0), 0);
	EXPECT_EQ(g.X(1,0), 0);
	EXPECT_EQ(g.X(2,0), 1);
	EXPECT_EQ(g.X(3,0), 1);
	EXPECT_EQ(g.X(0,1), 1);
	EXPECT_EQ(g.X(1,1), 1);
	EXPECT_EQ(g.X(2,1), 0);
	EXPECT_EQ(g.X(3,1), 0);
	// no basis por picking one or the other in the last SNP
	EXPECT_EQ(g.X(0,2), 1);
	EXPECT_EQ(g.X(1,2), 1);
	EXPECT_EQ(g.X(2,2), 0);
	EXPECT_EQ(g.X(3,2), 0);

	// co-dominant: the minor allele causes a phenotype when in heterozygosis
	CGWASDataHelper::encodeHeterozygousData(&g, 3);

	EXPECT_EQ(g.X(0,0), 0);
	EXPECT_EQ(g.X(1,0), 0);
	EXPECT_EQ(g.X(2,0), 1);
	EXPECT_EQ(g.X(3,0), 0);
	EXPECT_EQ(g.X(0,1), 0);
	EXPECT_EQ(g.X(1,1), 1);
	EXPECT_EQ(g.X(2,1), 0);
	EXPECT_EQ(g.X(3,1), 0);
	EXPECT_EQ(g.X(0,2), 0);
	EXPECT_EQ(g.X(1,2), 0);
	EXPECT_EQ(g.X(2,2), 0);
	EXPECT_EQ(g.X(3,2), 0);

}