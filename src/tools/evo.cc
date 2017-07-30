//
// Created by hclimente on 25/07/2017.
//

#include "gin/shake.h"
#include "boost/program_options.hpp"

int main(int argc, char* argv[]) {

	// settings
	namespace po = boost::program_options;
	po::options_description desc("Allowed options");

	desc.add_options()
			("ped,p", po::value<string>(), "Plink genotype file.")
			("pheno,f", po::value<string>(), "Plink phenotype file.")
			("net,n", po::value<string>(), "Sparse network file.")
			("association_score,t", po::value<string>()->default_value("chi2"), "Association score.")
			("model_selection,x", po::value<string>()->default_value("cons"), "Metric to evaluate models.")
			("depth,y", po::value<int>()->default_value(3), "Depth of the grid search.")
			("maf,m", po::value<double>()->default_value(0.05), "Float minor allele frequency filter.")
			("lambda,l", po::value<double>()->default_value(-1), "Lambda parameter.")
			("eta,e", po::value<double>()->default_value(-1), "Eta parameter.")
			("outdir,o", po::value<string>()->default_value("."), "Output directory.")
			("encoding,s", po::value<string>()->default_value("additive"), "snp_encoding.")
			("pc,c", po::value<int>()->default_value(0), "PC.")
			("seed,z", po::value<int>()->default_value(0), "Random state seed.")
			("debug,d", po::bool_switch()->default_value(false), "Debug flag (display extra information).")
			("help,h", "Produce this help message and exit.")
			;

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);

	if (vm.count("help")) {
		cout << desc << "\n";
		return 0;
	}

	double total_data = clock();
	double total = clock();
	string genotype_str = vm["ped"].as<string>();
	string phenotype_str = vm["pheno"].as<string>();
	string network_str = vm["net"].as<string>();
	double maf = vm["maf"].as<double>();
	string outfolder_str = vm["outdir"].as<string>();
	string snp_encoding = vm["encoding"].as<string>();
	int pcs = vm["pc"].as<int>();
	int seed = vm["seed"].as<int>();
	string association_score = vm["association_score"].as<string>();
	double lambda = vm["lambda"].as<double>();
	double eta = vm["eta"].as<double>();
	bool debug = vm["debug"].as<bool>();
	string model_selection = vm["model_selection"].as<string>();
	int depth = vm["depth"].as<int>();

	uint encoding = 0;
	if(snp_encoding == "additive") encoding = 0;
	else if(snp_encoding == "recessive") encoding = 1;
	else if(snp_encoding == "dominant") encoding = 2;
	else if(snp_encoding == "overdominant") encoding = 3;
	else {
		logging(ERROR,"Encoding does not exist!");
		exit(-1);
	}

	Shake shake = Shake();

	// read the data

	shake.readGWAS(genotype_str, encoding);
	shake.readNetwork(network_str);

	// Preprocess the data
	// shake.filterMAF(maf);
	// shake.adjustPC(pcs);

	VectorXd etas(4);
	etas << 0, 1, 2, 3;
	VectorXd lambdas(4);
	lambdas << 0, 1, 2, 3;

	shake.searchHyperparameters(10, model_selection, association_score);
	shake.selectSnps();

	shake.writeResults("results.txt");

}