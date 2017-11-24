#include "gin/gwas/CScones.h"
#include "gin/io/CSconesIO.h"
#include "math.h"
#include "gin/model_selection/CCrossValidation.h"
#include "gin/utils/CMatrixHelper.h"
#include "gin/utils/StringHelper.h"
#include "maxflow/maxflow.h"
/*
*Initialize default scones settings
*/
CSconesSettings::CSconesSettings() {
	folds = 10;
	seed = 0;
	selection_criterion = CONSISTENCY;
	selection_ratio = 0.8;
	test_statistic = SKAT;
	autoParameters = true;
	nParameters = 10;
	lambdas = VectorXd::Zero(nParameters);
	etas = VectorXd::Zero(nParameters);
	evaluateObjective = false;
	dump_intermediate_results = true;
	dump_path = "tmp/";
    gridsearch_depth = 1;
}

/*
* Default constructor for Scones
*/
CScones::CScones() {
	CSconesSettings settings;
	__settings = settings;
	__covs_set = false;
    __binary_y = false;
}

CScones::CScones(CSconesSettings const& settings) {
	__settings = settings;
	__covs_set = false;
    __binary_y = false;
}

CScones::CScones(VectorXd const& y, MatrixXd const& X, SparseMatrixXd const& L) throw (CSconesException) {
	CSconesSettings settings;
	__settings = settings;
	__covs_set = false;
    __binary_y = false;

	__y = y;
	__X = X;
	__W = L;

	__checkdata();

	__n_samples = X.rows();
	__n_features = X.cols();

	__sW = DiagXd(__n_features);
	__sW.diagonal() = VectorXd::Ones(__n_features);
	__L = __computeLaplacianMatrix();
}

CScones::CScones(VectorXd const& y, MatrixXd const& X, SparseMatrixXd const& L, CSconesSettings const& settings) throw (CSconesException) {
	__settings = settings;
	__covs_set = false;
    __binary_y = false;

	__y = y;
	__X = X;
	__W = L;

	__checkdata();

	__n_samples = X.rows();
	__n_features = X.cols();

	__sW = DiagXd(__n_features);
	__sW.diagonal() = VectorXd::Ones(__n_features);
	__L = __computeLaplacianMatrix();
}

CScones::CScones(VectorXd const& y, MatrixXd const& X, SparseMatrixXd const& L, MatrixXd const& covs) throw (CSconesException) {
	CSconesSettings settings;
	__settings = settings;
	__covs_set = true;
    __binary_y = false;

	__covs = covs;
	__y = y;
	__X = X;
	__W = L;

	__checkdata();

	__n_samples = X.rows();
	__n_features = X.cols();

	__sW = DiagXd(__n_features);
	__sW.diagonal() = VectorXd::Ones(__n_features);
	__L = __computeLaplacianMatrix();
}

CScones::CScones(VectorXd const& y, MatrixXd const& X, SparseMatrixXd const& L, MatrixXd const& covs, CSconesSettings const& settings) throw (CSconesException) {
	__settings = settings;
	__covs_set = true;
    __binary_y = false;

	__covs = covs;
	__y = y;
	__X = X;
	__W = L;

	__checkdata();

	__n_samples = X.rows();
	__n_features = X.cols();

	__sW = DiagXd(__n_features);
	__sW.diagonal() = VectorXd::Ones(__n_features);
	__L = __computeLaplacianMatrix();
}

void CScones::__checkdata() throw (CSconesException) {
	if(__y.cols()>1) throw CSconesException("Phenotype y has wrong dimensions! (n x 1)!");
	if(__X.rows() != __y.rows()) throw CSconesException("Genotype X and Phenotype y must have the same number of samples n!");
	if(__W.rows() != __W.cols()) throw CSconesException("Network L must be a squared matrix with size features m x m!");
	//if(__W.rows() != __X.cols()) throw CSconesException("Network L and Genotype X must have the same number of features m!");
	if(__covs_set==true) {
		if (__covs.rows()!=__X.rows()) throw CSconesException("Covariates and genotype must have the same number of samples n!");
	}
}

/*
 *Getter and Setter Methods
 */

//Set weight vector for SKAT association test
void CScones::setSKATWeights(VectorXd const& w) {
	if(w.rows()!=__X.cols()) throw CSconesException("Weight vector w has wrong dimensions");
	__sW = DiagXd(__n_features);
	__sW.diagonal() = VectorXd::Ones(__n_features).cwiseProduct(w);
}

VectorXd CScones::getIndicatorVector() {
	return __indicator_vector;
}

float64 CScones::getBestLambda() {
	return __best_lambda;
}

float64 CScones::getBestEta() {
	return __best_eta;
}

SparseMatrixXd CScones::getW() {
    return __W;
}

MatrixXd CScones::getCMatrix() {
	return __cMat;
}
CSconesSettings CScones::getSettings() {
    return __settings;
}

std::vector<std::vector<std::vector<SparseMatrixXd> > > CScones::getResultStack() {
	return __result_stack;
}

float64 CScones::getObjectiveScore() {
    return __objective_score;
}

/*
* Main algorithm
*/
void CScones::__selectRegressionModel() {
	bool binary = true;
	for(int64 i=0; i<__y.rows();i++) {
		if(!(__y(i)==0 || __y(i)==1)) {
			binary = false;
			break;
		}
	}
	if(binary) { //if phenotype is binary select LogisticRegression
        __logistic_regression = CLogisticRegression();
        __binary_y = true;
	} else { //if phenotype is continuous select LinearRegression
		__linear_regression = CLinearRegression();
        __binary_y = false;
	}
}

VectorXd CScones::__computeSKATScore(MatrixXd const& X, VectorXd const& r) {
	VectorXd nonweighted_skat = ((X.transpose()*r).array().pow(2));
	return __sW*nonweighted_skat;
}

VectorXd CScones::__computeChisqScore(MatrixXd const& X, VectorXd const& r) {
    VectorXd chisq(X.cols());

    for (int i = 0; i < X.cols(); i++) {
        MatrixXd tab = CChi2::get2DContingencyTable(r, X.col(i), true);
        chisq(i) = CChi2::calculateChi2(tab);
    }

    return chisq;

}

VectorXd CScones::__computeCochranArmitageT(MatrixXd const& X, VectorXd const& r, string const& geneticModel) {
    VectorXd T(X.cols());
    VectorXd model(3);

    if(geneticModel == "dominance")
		model << 1, 1, 0;
	else if(geneticModel == "recessive")
		model << 0, 1, 1;
	else if(geneticModel == "codominance")
		model << 0, 1, 2;
	else
        model << 1, 1, 1;

    for (int i = 0; i < X.cols(); i++) {
		MatrixXd tab = CChi2::get2DContingencyTable(r, X.col(i), true);
        T(i) = CChi2::calculateChi2Trend(tab, model);
    }

    return T;
}

VectorXd CScones::__computeScoreStatistic(MatrixXd const& X, VectorXd const& r) {
	VectorXd score;
	if(__settings.test_statistic==SKAT)
        score = __computeSKATScore(X,r);
    else if(__settings.test_statistic==CHISQ)
        score = __computeChisqScore(X,r);
    else if(__settings.test_statistic==TREND)
        score = __computeCochranArmitageT(X, r, "dominance");
	return score;
}

void CScones::__autoParameters() {
	VectorXd r;
	//fit linear regression without evaluating the likelihood function
	if(__covs_set) {
		if(__binary_y==true) {
            __logistic_regression.fit(__y,__covs);
            r = __logistic_regression.getResiduals();
        } else {
            __linear_regression.fit(__y,__covs);
            r = __linear_regression.getResiduals();
	    }
    } else { //if no covariates are set use the zero meaned y ase residuals
		r = __y.array() - __y.mean();
	}
	//set parameters vectors
    VectorXd c = __computeScoreStatistic(__X,r);
	VectorXd tmpc = (c.array()==0).select(c.maxCoeff(),c);
	float64 minc = floor(log10(tmpc.minCoeff()));
	float64 maxc = ceil(log10(c.maxCoeff()));

	__settings.etas = VectorXd::LinSpaced(__settings.nParameters, minc, maxc);
	for(int i=0; i<__settings.etas.rows(); i++)
		__settings.etas(i) = pow(10,__settings.etas(i));

	__settings.lambdas = VectorXd::LinSpaced(__settings.nParameters, minc - 1, maxc + 1);
	for(int i=0; i<__settings.lambdas.rows(); i++)
		__settings.lambdas(i) = pow(10,__settings.lambdas(i));
}

//void CScones::maxflow(SparseMatrixXd const& A, MatrixXd const& T, VectorXd* indicator_vector) {
void CScones::maxflow(SparseMatrixXd const &A, MatrixXd const &T, VectorXd *indicator_vector) {
	indicator_vector->resize(__n_features);
	//create graph out of adjacency matrix
	typedef Graph<float64, float64, float64> MaxGraph;
	MaxGraph *g = new MaxGraph(A.rows(),A.nonZeros());
	//Initialize nodes
	g->add_node(indicator_vector->rows());
	//traverse the sparse adjacency matrix A
	for(int64 k=0; k<A.outerSize(); k++) {
		for(Eigen::SparseMatrix<float64>::InnerIterator it(A,k); it; ++it) {
			g->add_edge(it.row(),it.col(),it.value(),0.0);
		}
	}
	//traverse the T matrix
	for(int64 i=0;i<2;i++) {
		for(int64 k=0; k<T.rows(); k++) {
			if(i==0) {
				if(T(k,i)!=0) {
					g->add_tweights(k,T(k,i),0.0);
				}
			} else {
				if(T(k,i)!=0) {
					g->add_tweights(k,0.0,T(k,i));
				}
			}
		}
	}
	//Run maxflow algorithm
	g->maxflow();

	//create indicator_vector
	for(int64 i=0; i<indicator_vector->rows();i++)
		(*indicator_vector)(i) = g->what_segment(i);
	//delete graph
	delete g;
}

void CScones::__optimize_objective(VectorXd const& c, float64 const& lambda, VectorXd* indicator_vector, float64* objective_score) {
	// Main Graph W: edges multiplied by lambda
	SparseMatrixXd lW = lambda * __W;
	//Add source and sink
	//if((c.array()>0).count()==0 || (c.array()<=0).count()==0)
	//	logging(WARNING,"WARNING: Trivial solution ahead!");

	// Matrix A containing As and At
	MatrixXd A(__n_features,2);
	//connect positive c values to sink
	VectorXd pos_c = (c.array() <= 0).select(0, c);
	//connect negative c values to source
	VectorXd neg_c = -c;
	neg_c = (c.array() > 0).select(0, neg_c);
	//Store data
	A.col(0) = neg_c;
	A.col(1) = pos_c;

	//compute maxflow
	maxflow(lW, A, indicator_vector);

	//compute objective score
	if(__settings.evaluateObjective) {
		(*objective_score) = 0.0;
		for(uint i=0; i<c.rows(); i++) {
			if((*indicator_vector)(i)==0) {
				(*objective_score) += c(i);
			}
		}

		for(int64 k=0; k<lW.outerSize(); k++) {
			for(Eigen::SparseMatrix<float64>::InnerIterator it(lW,k); it; ++it) {
				if((*indicator_vector)(it.row())==1)
					if((*indicator_vector)(it.col())==0)
						(*objective_score) += it.value();
			}
		}
		(*objective_score) -= c.sum();
	} else {
		(*objective_score) = NAN;
	}
}

void CScones::__gridsearch(VectorXd const& y, MatrixXd const& x, MatrixXd const& covs) throw (CSconesException) {
	VectorXd r;
	if(__covs_set) {
		if(__binary_y==true) {
            __logistic_regression.fit(y,covs);
            r = __logistic_regression.getResiduals();
        } else {
            __linear_regression.fit(y,covs);
            r = __linear_regression.getResiduals();
	    }
	} else {
		r = y.array() - y.mean();
	}
	VectorXd c = __computeScoreStatistic(x,r);
	VectorXd c_transformed;
	std::vector<std::vector<SparseMatrixXd> > lambda_stack;
	//Perform gridsearch to find optimal lambda and eta values
	for(uint i=0; i<__settings.lambdas.rows(); i++) {
		std::vector<SparseMatrixXd> eta_stack;
		for(uint j=0;j<__settings.etas.rows(); j++) {
			//Substract eta value from score-statistic c
			c_transformed = c.array()-__settings.etas(j);
			//Optimize objective
			VectorXd indicator_vector;
			float64 objective_score;
			__optimize_objective(c_transformed, __settings.lambdas(i), &indicator_vector, &objective_score);
			eta_stack.push_back(indicator_vector.sparseView());
		}
		lambda_stack.push_back(eta_stack);
	}
	__result_stack.push_back(lambda_stack);
}

void CScones::test_associations() throw (CSconesException) {
	//Select the correct Regression Model
	__selectRegressionModel();
	if(__settings.autoParameters) {
		__autoParameters();
	} else {
		if(__settings.etas.rows() == 0 || __settings.lambdas.rows() == 0)
			throw CSconesException("Array of lambda or eta values cannot be empty!");
	}

    //Perform crossvalidated gridsearch to find parameters
    CCrossValidation cv(__settings.seed);
    cv.kFold(__settings.folds,__n_samples);
    MatrixXd::Index best_eta_index, best_lambda_index;

	logging(STATUS, "Choosing eta and lambda values.");

    for (int i=1; i <= __settings.gridsearch_depth; i++){

		logging(INFO, "Grid search " + StringHelper::to_string<int>(i) + " / " + StringHelper::to_string<int>(__settings.gridsearch_depth) + ".");

        for(uint k=0;k<__settings.folds;k++) {
            VectorXd tr_indices = cv.getTrainingIndices(k);
            //Slice data according to the indices
            MatrixXd x_train = sliceRowsMatrix(__X,tr_indices);
            VectorXd y_train = sliceRowsMatrix(__y,tr_indices);
            MatrixXd cov_train;
            if(__covs_set) {
                cov_train = sliceRowsMatrix(__covs,tr_indices);
            }
            __gridsearch(y_train,x_train,cov_train);
        }

        //Evaluate all solutions and select the best eta and lambda according to the selection criterion
        if(__settings.selection_criterion==CONSISTENCY)
        {
            __cMat = __evaluateConsistency();
            __best_c = __cMat.maxCoeff(&best_eta_index,&best_lambda_index);
        }
        else if(__settings.selection_criterion == AICc | __settings.selection_criterion == BIC | __settings.selection_criterion == mBIC)
        {
            __cMat = __evaluateInformation();
            __best_c = __cMat.minCoeff(&best_eta_index,&best_lambda_index);
        }

        __best_lambda = __settings.lambdas[best_lambda_index];
        __best_eta = __settings.etas[best_eta_index];

        // if it's not the last iteration,
        // recalculate lambda and eta ranges according to best lambda and eta found
        if (i < __settings.gridsearch_depth){
            float64 order_lambdas = log10(__best_lambda);
            float64 incr_lambdas = (log10(__settings.lambdas.maxCoeff()) - log10(__settings.lambdas.minCoeff())) * 0.2;
            __settings.lambdas = VectorXd::LinSpaced(__settings.nParameters, order_lambdas - incr_lambdas, order_lambdas + incr_lambdas);
            for(int i=0; i<__settings.lambdas.rows(); i++)
                __settings.lambdas(i) = pow(10,__settings.lambdas(i));

            float64 order_etas = log10(__best_eta);
            float64 incr_etas = (log10(__settings.etas.maxCoeff()) - log10(__settings.etas.minCoeff())) * 0.2;
            __settings.etas = VectorXd::LinSpaced(__settings.nParameters, order_etas - incr_etas, order_etas + incr_etas);
            for(int i=0; i<__settings.etas.rows(); i++)
                __settings.etas(i) = pow(10,__settings.etas(i));

            // clear __result_stack for next iteration
            __result_stack.clear();
        }
    }

	logging(STATUS, "Feature selection.");

    //Selected Features for the best eta and lambda are those selected in all folds
    __indicator_vector = VectorXd::Zero(__n_features);
    if(__best_lambda>0.0 && __best_eta>0) {
        for(uint k=0; k<__settings.folds; k++) {
            __indicator_vector += __result_stack[k][best_lambda_index][best_eta_index];
        }
        //Select only those features selected in all folds
        __indicator_vector /= __settings.folds;
        for(uint64 i=0; i<__indicator_vector.rows();i++) {
            if(__indicator_vector(i)!=1) __indicator_vector(i)=0;
        }
    }

    /*ALTERNATIVE TO GET THE FINAL INDICATOR VECTOR, retrain model using the optimized parameters:
     *     test_associations(__best_lambda,__best_eta);
     * HOWEVER this sometimes leads to different selections when some solutions are not picked
     * in all folds but are in the global training. */
}

MatrixXd CScones::__evaluateConsistency() throw (CSconesException){
    MatrixXd gridResults;
    // N number of features
    float64 N = __W.outerSize();
    float64 maxn = ceil(N*0.01);
    gridResults = MatrixXd::Zero(__settings.etas.rows(),__settings.lambdas.rows());
    for(int e=0; e<__settings.etas.rows();e++){
        for(int l=0; l<__settings.lambdas.rows();l++) {
            float64 cindex = 0.0;
            for(uint k1=0; k1<__settings.folds; k1++) {
                //inds.push_back(__result_stack[k1][l][e]);
                // count number of non-zero coefficients in sparse matrix
                float64 n_k1 = __result_stack[k1][l][e].nonZeros();
                // if trivial solutions (no SNP/all SNPs/more than a threshold), skip this fold
                if (n_k1==0 || n_k1 == N || n_k1 > maxn) continue;
                // else scan what happens in folds with an index > than the current one
                for(uint k2=k1+1; k2<__settings.folds; k2++) {
                    float64 n_k2 = __result_stack[k2][l][e].nonZeros();
                    if (n_k2==0 || n_k2 == N || n_k2 > maxn) continue;

                    // normalize the consistency using the maximum possible consistency ie max 1
                    float64 consistency = N * (__result_stack[k1][l][e].cwiseProduct(__result_stack[k2][l][e])).sum() - n_k1 * n_k2;
                    float64 maxConsistency = N * fmin(n_k1, n_k2) - n_k1 * n_k2;
                    cindex += consistency/maxConsistency;
                }
            }
            cindex = 2.0 * cindex / (__settings.folds * (__settings.folds - 1));
            gridResults(e,l) = cindex;
        }
    }

    return gridResults;
}

MatrixXd CScones::__evaluateInformation() throw (CSconesException) {
    MatrixXd gridResults;

    // number of features
    float64 N = __W.outerSize();
    float64 maxn = ceil(N*0.1);
    gridResults = MatrixXd::Zero(__settings.etas.rows(),__settings.lambdas.rows());
    for(int e=0; e<__settings.etas.rows();e++){

        for(int l=0; l<__settings.lambdas.rows();l++) {
            VectorXd indicator_vector = VectorXd::Zero(__n_features);
            for(uint k = 0; k < __settings.folds; k++){
                float64 n = __result_stack[k][l][e].nonZeros();
                // if trivial solutions (no SNP/all SNPs/more than a threshold), skip this fold
                if (n==0 || n == N || n > maxn) continue;
                indicator_vector += __result_stack[k][l][e];
            }

            // take only features picked across all folds
            indicator_vector /= __settings.folds;
            for(uint64 i=0; i<indicator_vector.rows();i++) {
                if(indicator_vector(i)!=1) indicator_vector(i)=0;
            }

            // fit a model to the selected features
            // max possible value if no features selected
            float64 informationMetric = 1e31;
            
            if (indicator_vector.sum() != 0)
			{
                MatrixXd x_tr = sliceColsMatrixByBinaryVector(__X, indicator_vector);
                if(__binary_y==true) {
                    __logistic_regression.fit(__y, x_tr);
					if(__settings.selection_criterion == AICc)
						informationMetric = __logistic_regression.getAICc();
					else if(__settings.selection_criterion == BIC)
						informationMetric = __logistic_regression.getBIC();
					else if(__settings.selection_criterion == mBIC)
						informationMetric = __logistic_regression.getAICc() - indicator_vector.transpose() * __L * indicator_vector;
                } else {
                    __linear_regression.fit(__y, x_tr);
					if(__settings.selection_criterion == AICc)
						informationMetric = __linear_regression.getAICc();
					else if(__settings.selection_criterion == BIC)
						informationMetric = __linear_regression.getBIC();
					else if(__settings.selection_criterion == mBIC)
						informationMetric = __logistic_regression.getAICc() - indicator_vector.transpose() * __L * indicator_vector;
                }
            }
            gridResults(e,l) = informationMetric;
        }
    }

    return gridResults;
}

void CScones::test_associations(float64 const& lambda, float64 const& eta) {
	//Select the correct Regression Model
	__selectRegressionModel();
	//Test associations for a fixed lambda and eta
	VectorXd r;
	if(__covs_set) {
		if(__binary_y==true) {
            __logistic_regression.fit(__y,__covs);
            r = __logistic_regression.getResiduals();
        } else {
            __linear_regression.fit(__y,__covs);
            r = __linear_regression.getResiduals();
	    }
	} else {
		r = __y.array() - __y.mean();
	}
	VectorXd c = __computeScoreStatistic(__X,r).array() - eta;
	//solve objective function
	__optimize_objective(c, lambda, &__indicator_vector, &__objective_score);
}

VectorXd CScones::getObjectiveFunctionTerms(float64 const& lambda, float64 const& eta){
	VectorXd terms(3);

	double connectivity = lambda * __indicator_vector.transpose() * __L * __indicator_vector;
	double sparsity = eta * __indicator_vector.sum();

    // association
    double association = 0;
	VectorXd c = getScoreStatistic();

	for(uint i=0; i<__indicator_vector.size(); i++) {
		if(__indicator_vector(i)!=0) {
			association += c(i);
		}
	}

	terms << association, connectivity, sparsity;

	return terms;
}

VectorXd CScones::getScoreStatistic() {
	VectorXd r;
	if(__covs_set) {
		if(__binary_y) {
			r = __logistic_regression.getResiduals();
		} else {
			r = __linear_regression.getResiduals();
		}
	} else {
		r = __y.array() - __y.mean();
	}

	VectorXd c = __computeScoreStatistic(__X,r);

	return c;
}

SparseMatrixXd CScones::__computeLaplacianMatrix() {

	SparseMatrixXd W = __W;
	SparseMatrixXd D(__n_features, __n_features);

	typedef Eigen::Triplet<double> T;
	std::vector<T> diagonal;
	diagonal.reserve(__n_features);

	VectorXd degree = VectorXd::Zero(__n_features);
	for (int k = 0; k < W.outerSize(); ++k)
		for (SparseMatrixXd::InnerIterator it(W, k); it; ++it)
			degree(it.row()) += 1;

	for (int i = 0; i < __n_features; i++)
		diagonal.push_back(T(i, i, degree(i)));

	D.setFromTriplets(diagonal.begin(), diagonal.end());
	SparseMatrixXd L = D;
	L -= W;

	return L;
}
