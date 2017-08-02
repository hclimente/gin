//
// Created by hclimente on 21/07/2017.
//

#ifndef GIN_SHAKE_H
#define GIN_SHAKE_H

#include "gin/model_selection/grid_cv.h"
#include "gin/feature_selection/scones.h"
#include "gin/stats/univariate_association.h"

#include "gin/gwas/CGWASData.h"
#include "gin/io/CPlinkParser.h"
#include "gin/io/CSconesIO.h"

class Shake {
public:
	Shake();
	~Shake();

	// GWAS operations
	void readGWAS(string const&, uint);
	void readNetwork(string const&);
	void searchHyperparameters(uint, uint const&, uint const&);
	void selectSnps();
	void writeResults(string const&);

	// setters & getters
	VectorXd selectedSnps() { return __selectedSnps; }
	double bestEta() { return __bestEta; }
	double bestLambda() { return __bestLambda; }
	GWASData* gwas() { return __gwas; }

private:
	GWASData* __gwas;
	GridCV* __cvgrid;

	double __bestEta;
	double __bestLambda;
	VectorXd __c;
	VectorXd __selectedSnps;

};

#endif //GIN_SHAKE_H
