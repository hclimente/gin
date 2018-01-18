//
// Created by hclimente on 25/07/2017.
//

#ifndef GIN_SETTINGS_H
#define GIN_SETTINGS_H

#include "gin/globals.h"

#ifndef AS_GINLIB

#include "boost/program_options.hpp"
namespace po = boost::program_options;

#endif

class Settings
{
public:

	Settings(int argc, const char* const argv[]) {
		__error = false;
		__parseOptions(argc, argv);
	}

	Settings(std::string pedBasename, std::string networkFilename, uint encoding, uint modelScore, uint associationScore, VectorXd eta, VectorXd lambda, std::string output) {
		__error = false;
		__debug = false;
		__pedBasename = pedBasename;
		__networkFilename = networkFilename;
		__encoding = encoding;
		__modelScore = modelScore;
		__associationScore = associationScore;
		__etas = eta;
		__lambdas = lambda;
		__output = output;
	}

	void printHelp() {
		std::cout << __options() << "\n";
	}

	std::string pedBasename() { return __pedBasename; }
	uint encoding() { return __encoding; }
	std::string networkFilename() { return __networkFilename; }
	uint modelScore() { return __modelScore; }
	uint associationScore() { return __associationScore; }
	VectorXd etas() { return __etas; }
	VectorXd lambdas() { return __lambdas; }
	std::string output() { return __output; }
	bool debug() { return __debug; }
	bool error() { return __error; }

private:

	std::string __pedBasename;
	std::string __networkFilename;
	// double __maf;
	std::string __output;
	// int __pcs;
	// int __seed;
	uint __associationScore;
	VectorXd __lambdas;
	VectorXd __etas;
	bool __debug;
	uint __modelScore;
	// int __depth;
	uint __encoding;

	bool __error;

	#ifndef AS_GINLIB

	po::options_description __options() {
		po::options_description desc("Allowed options");

		desc.add_options()
				("ped,p", po::value<std::string>(), "Plink genotype file.")
				("net,n", po::value<std::string>(), "Sparse network file.")
				("association_score,c", po::value<std::string>()->default_value("chi2"), "Association score.")
				("model_score,m", po::value<std::string>()->default_value("cons"), "Metric to evaluate the models.")
				// ("depth,y", po::value<int>()->default_value(3), "Depth of the grid search.")
				// ("maf,m", po::value<double>()->default_value(0.05), "Float minor allele frequency filter.")
				("lambda,l", po::value<double>()->default_value(-1), "Lambda parameter.")
				("eta,e", po::value<double>()->default_value(-1), "Eta parameter.")
				("out,o", po::value<std::string>()->default_value("results.txt"), "Output file.")
				("encoding,e", po::value<std::string>()->default_value("additive"), "Model of inheritance (additive, recessive, dominant or codominant).")
				// ("pc,c", po::value<int>()->default_value(0), "PC.")
				// ("seed,z", po::value<int>()->default_value(0), "Random state seed.")
				("debug,d", po::bool_switch()->default_value(false), "Debug flag (display extra information).")
				("help,h", "Produce this help message and exit.")
				;

		return desc;
	}

	void __parseOptions(int argc, const char* const argv[]) {

		po::variables_map vm;
		po::store(po::parse_command_line(argc, argv, __options()), vm);
		po::notify(vm);

		if (vm.count("help")) {
			__error = true;
		}

		__pedBasename = vm["ped"].as<std::string>();
		__networkFilename = vm["net"].as<std::string>();
		// __maf = vm["maf"].as<double>();
		__output = vm["out"].as<std::string>();
		std::string snp_encoding = vm["encoding"].as<std::string>();
		// __pcs = vm["pc"].as<int>();
		// __seed = vm["seed"].as<int>();
		std::string association_score_str = vm["association_score"].as<std::string>();
		double lambda = vm["lambda"].as<double>();
		double eta = vm["eta"].as<double>();
		__debug = vm["debug"].as<bool>();
		std::string model_selection_str = vm["model_score"].as<std::string>();
		// __depth = vm["depth"].as<int>();

		__encoding = 0;
		if(snp_encoding == "additive") {
			__encoding = 0;
		} else if(snp_encoding == "recessive") {
			__encoding = 1;
		} else if(snp_encoding == "dominant") {
			__encoding = 2;
		} else if(snp_encoding == "codominant") {
			__encoding = 3;
		} else {
			logging(ERROR,"Encoding does not exist!");
			__error = true;
		}

		__associationScore = SKAT;
		if (association_score_str == "chi2") {
			__associationScore = CHI2;
		} else if (association_score_str == "trend") {
			__associationScore = TREND;
		}

		__modelScore = CONSISTENCY;
		if (model_selection_str == "bic") {
			__modelScore = BIC;
		} else if (model_selection_str == "aic") {
			__modelScore = AIC;
		} else if (model_selection_str == "aicc") {
			__modelScore = AICc;
		} else if (model_selection_str == "mbic") {
			__modelScore = mBIC;
		}

		if (eta != -1 & lambda != -1) {
			__etas = VectorXd(1);
			__etas(0) = eta;
			__lambdas = VectorXd(1);
			__lambdas(0) = lambda;
		}


	}

	#endif //AS_GINLIB

};

#endif //GIN_SETTINGS_H
