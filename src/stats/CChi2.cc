#include "gin/stats/CChi2.h"
//#include "CGamma.h"

void CChi2::__checkParameters(float64 const& k) throw (CChi2Exception) {
	if (k==0) throw CChi2Exception("Degress of freedom (k) cannont be zero!");
	else if (k<0) throw CChi2Exception("Degress of freedom (k) cannont be negative!");
}

float64 CChi2::cdf(float64 const& x, float64 const& k) throw (CChi2Exception) {
	__checkParameters(k);
    return Cephes::chdtr(k,x);
    //Own Implementation
    /*if (k==2.0) { //Special case if k==2.0
		return 1.0 - exp(-0.5*x);
	} else { //for any k
		return (CGamma::Special::regularizedLowerIncompleteGamma(0.5*x,0.5*k));
	}*/
}

float64 CChi2::logcdf(float64 const& x, float64 const& k) throw (CChi2Exception) {
	__checkParameters(k);
	return log(cdf(x,k));
}

float64 CChi2::sf(float64 const& x, float64 const& k) throw (CChi2Exception) {
	__checkParameters(k);
    if(x < 0) return 1;
    else return Cephes::chdtrc(k,x);
    //return 1-cdf(x,k); //is too less accurate! Use complemented approximation
    //OWN IMPlEMENTATION
	//return CGamma::Special::complementedIncompleteGamma(0.5*x,0.5*k);
}

float64 CChi2::isf(float64 const& x, float64 const& k) throw (CChi2Exception) {
	__checkParameters(k);
	if(x>1.0 || x<0.0) throw CChi2Exception("x has to be in the range [0.0,1.0]!");
    return Cephes::chdtri(k,x);
}

float64 CChi2::logsf(float64 const& x, float64 const& k) throw (CChi2Exception) {
	__checkParameters(k);
	return log(sf(x,k));
}

float64 CChi2::pdf(float64 const& x, float64 const& k) throw (CChi2Exception) {
	__checkParameters(k);
	if(x<0.0) return 0.0;
	return pow(x,0.5*k-1)*exp(-0.5*x)/(pow(2.0,0.5*k)*tgamma(0.5*k));
}

float64 CChi2::logpdf(float64 const& x, float64 const& k) throw (CChi2Exception) {
	__checkParameters(k);
	return log(pdf(x,k));
}

MatrixXd CChi2::get2DContingencyTable(VectorXd const& x, VectorXd const& Y, bool pseudoCounts) throw (CChi2Exception) {

	if (x.size() != Y.size()) throw CChi2Exception("Variable vector lengths are different.");

	std::map<int, std::map<int,int> > counts;

	for (int i = 0; i < x.size(); i++)
		counts[Y(i)][x(i)] += 1;

	MatrixXd table = MatrixXd::Zero(2, 3);
	if (pseudoCounts) {
		table = MatrixXd::Ones(2, 3);
	}

	int row = 0;
	for (std::map<int, std::map<int,int> >::iterator it = counts.begin(); it!= counts.end(); ++it) {
		int col = 0;
		for (int i = 0; i <= 2; i++){

			std::map<int,int>::iterator f = it->second.find(i);
			if (f != it->second.end()) {
				table(row, col) = table(row, col) + it->second[i];
			}
			col++;
		}
		row++;
	}

	return table;
}

double CChi2::calculateChi2(MatrixXd const& table) throw (CChi2Exception) {
	double N = table.sum();

	MatrixXd p_phenotype = table.rowwise().sum() / N;
	MatrixXd p_genotype = table.colwise().sum() / N;
	MatrixXd p(table.rows(), table.cols());
	p = p_phenotype * p_genotype;

	MatrixXd numerator = (table/N - p).array().square();

	float64 chisq = N * numerator.cwiseQuotient(p).sum();

	return chisq;
}

float64 CChi2::calculateChi2Trend(MatrixXd const& table, VectorXd const& model) throw (CChi2Exception) {
	double N = table.sum();

	double T1 = (table.row(0) * model).sum();
	double T2 = (table.colwise().sum() * model).sum();
	double T3 = (table.colwise().sum() * model.array().pow(2).matrix()).sum();
	long int C = table.row(0).sum();
	long int D = table.row(1).sum();
	double V = C * D * (N * T3 - pow(T2, 2)) / (pow(N, 2) * (N - 1));

	double chisq = pow(T1 - C * T2 / N, 2) / V ;

	return chisq;
}
