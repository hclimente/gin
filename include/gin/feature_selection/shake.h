//
// Created by hclimente on 21/07/2017.
//

#ifndef GIN_SHAKE_H
#define GIN_SHAKE_H

#include "gin/io/grid_views.h"
#include "gin/model_selection/grid_cv.h"
#include "gin/feature_selection/scones.h"
#include "gin/stats/univariate_association.h"

#include "gin/gwas/CGWASData.h"
#include "gin/io/CPlinkParser.h"
#include "gin/io/CSconesIO.h"
#include "gin/utils/StringHelper.h"

class Shake {
public:
	Shake();
	Shake(MatrixXd*, VectorXd*, SparseMatrixXd*);
	~Shake();

	// GWAS operations
	void readGWAS(string const&, uint);
	void readNetwork(string const&);
	void selectHyperparameters(uint, uint const &, uint const &);
	void selectHyperparameters(uint, uint const &, uint const &, VectorXd const &, VectorXd const&);
	void selectSNPs();
	void writeResults(string const&);

	// setters & getters
	VectorXd selectedSnps() { return __selectedSnps; }
	double bestEta() { return __bestEta; }
	double bestLambda() { return __bestLambda; }
	VectorXd c() { return __c; }
	GWASData* gwas() { return __gwas; }
	GridCV* grid() { return __cvgrid; }

	void setDebug(bool debug) { __debug = debug; }

private:
	GWASData* __gwas;
	GridCV* __cvgrid;

	double __bestEta;
	double __bestLambda;
	VectorXd __c;
	VectorXd __selectedSnps;
	bool __debug;

	void __computeAssociation(uint const &);

};

#endif //GIN_SHAKE_H
