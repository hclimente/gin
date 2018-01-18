#include <iostream>
#include "boost/program_options.hpp"

#include "gin/gwas/CScones.h"
#include "gin/io/CPlinkParser.h"
#include "gin/utils/CKernels.h"
#include "gin/io/CSconesIO.h"
#include "gin/stats/CStats.h"
#include "gin/utils/CMatrixHelper.h"

using namespace std;

int main(int argc, char* argv[]) {
	//get command line arguments
	// Declare the supported options.
	namespace po = boost::program_options;
	po::options_description desc("Allowed options");

	desc.add_options()
		("ped,p", po::value<string>(), "Plink genotype file.")
		("pheno,f", po::value<string>(), "Plink phenotype file.")
		("net,n", po::value<string>(), "Sparse network file.")
		("association_score,t", po::value<string>()->default_value("skat"), "Association score.")
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
	if(snp_encoding=="additive") encoding = 0;
	else if(snp_encoding=="recessive") encoding = 1;
	else if(snp_encoding=="dominant") encoding = 2;
	else if(snp_encoding=="overdominant") encoding = 3;
	else {
		logging(GIN_ERROR,"Encoding does not exist!");
		exit(-1);
	}

    CSconesSettings settings;
	settings.seed = seed;
    if (lambda != -1 & eta != -1){
        VectorXd l(1);
        l(0) = lambda;
        VectorXd e(1);
        e(0) = eta;
        settings.lambdas = l;
        settings.etas = e;
        // avoid gridsearch
        settings.autoParameters = false;
    }

    if (association_score == "chisq")
        settings.test_statistic = CHISQ;

	settings.gridsearch_depth = depth;

	if (model_selection == "cons") settings.selection_criterion = CONSISTENCY;
	else if (model_selection == "aicc") settings.selection_criterion = AICc;
	else if (model_selection == "bic") settings.selection_criterion = BIC;
	else if (model_selection == "mbic") settings.selection_criterion = mBIC;

	GWASData data;

	float64 begin;
	begin = clock();
	logging(GIN_STATUS,"Reading Genotype file...");
	CPlinkParser::readPEDFile(genotype_str + ".ped", &data);
	logging(GIN_INFO,"Number of SNPs: " + StringHelper::to_string<uint64>(data.n_snps));
	logging(GIN_INFO,"Number of Samples: " + StringHelper::to_string<uint64>(data.n_samples));
	logging(GIN_WARNING,"Finished in " + StringHelper::to_string<float64>(float64(clock()-begin)/CLOCKS_PER_SEC) + " sec\n");

	begin = clock();
	logging(GIN_STATUS,"Reading Mapping file...");
	CPlinkParser::readMAPFile(genotype_str + ".map", &data);
	logging(GIN_WARNING,"Finished in " + StringHelper::to_string<float64>(float64(clock()-begin)/CLOCKS_PER_SEC) + " sec\n");


	begin = clock();
	logging(GIN_STATUS,"Reading Phenotype file...");
	CPlinkParser::readPhenotypeFile(phenotype_str,&data);
	logging(GIN_WARNING,"Finished in " + StringHelper::to_string<float64>(float64(clock()-begin)/CLOCKS_PER_SEC) + " sec\n");


	begin = clock();
	logging(GIN_STATUS,"Encoding SNP data...");
	CGWASDataHelper::encodeHeterozygousData(&data,encoding);
	logging(GIN_WARNING,"Finished in " + StringHelper::to_string<float64>(float64(clock()-begin)/CLOCKS_PER_SEC) + " sec\n");

	/*
	begin = clock();
	logging(GIN_STATUS,"Filter unique SNPs...");
	CGWASDataHelper::filterUniqueSNPs(&data);
	logging(GIN_WARNING,"Finished in " + StringHelper::to_string<float64>((clock()-begin)/CLOCKS_PER_SEC) + " sec\n");
	*/

	begin = clock();
	logging(GIN_STATUS,"Filter SNPs by MAF...");
	CGWASDataHelper::filterSNPsByMAF(&data,maf);
	logging(GIN_INFO,"Number of SNPs (after MAF): " + StringHelper::to_string<uint64>(data.n_snps));
	logging(GIN_WARNING,"Finished in " + StringHelper::to_string<float64>(float64(clock()-begin)/CLOCKS_PER_SEC) + " sec\n");

	begin = clock();
	logging(GIN_STATUS,"Loading and filtering network file...");
	CSconesIO::readSparseNetworkFile(network_str,&data);
	logging(GIN_WARNING,"Finished in " + StringHelper::to_string<float64>(float64(clock()-begin)/CLOCKS_PER_SEC) + " sec\n");


	for(uint i=0; i<data.phenotype_names.size(); i++) {

		begin = clock();
		logging(GIN_STATUS,"Remove samples with missing values...");
		GWASData tmpData = CGWASDataHelper::removeSamples4MissingData(data,i);
		if(tmpData.n_samples!=data.n_samples) {
			logging(GIN_INFO,"#Samples removed: " + StringHelper::to_string<uint64>(data.n_samples-tmpData.n_samples));
		}
		logging(GIN_WARNING,"Finished in " + StringHelper::to_string<float64>(float64(clock()-begin)/CLOCKS_PER_SEC) + " sec\n");
		logging(GIN_WARNING,"Total Data: " + StringHelper::to_string<float64>(float64(clock()-total_data)/CLOCKS_PER_SEC) + " sec\n");

		//OPTIONAL SET DIFFERENT SConES Settings
		begin = clock();
		//COMPUTE PRINCIPLE COMPONENTS
		if (pcs>0) {
			begin = clock();
			logging(GIN_STATUS,"Computing Realized Relationship Kernel and Principle Components...");
			data.K = CKernels::realizedRelationshipKernel(data.X);
			MatrixXd PCs = CStats::principle_components(data.K);
			PCs = sliceColsMatrix(PCs,VectorXd::LinSpaced(pcs,0,pcs-1));
			logging(GIN_WARNING,"Finished in " + StringHelper::to_string<clock_t>((clock()-begin)/CLOCKS_PER_SEC) + " sec\n");

			logging(GIN_STATUS,"Running Scones for phenotype: " + data.phenotype_names[i]);
            CScones scones(tmpData.Y.col(i),tmpData.X,tmpData.network,PCs,settings);
            scones.test_associations();
			logging(GIN_WARNING,"Finished in " + StringHelper::to_string<float64>(float64(clock()-begin)/CLOCKS_PER_SEC) + " sec\n");

			begin = clock();
			string output_str = outfolder_str + "/" + data.phenotype_names[i] + ".scones.out.txt";
			logging(GIN_STATUS,"Writing output to " + output_str);
            CSconesIO::writeOutput(output_str, tmpData, scones.getIndicatorVector(),scones.getBestLambda(),scones.getBestEta());
			output_str = outfolder_str + "/" + data.phenotype_names[i] + ".scones.pmatrix.txt";
			logging(GIN_STATUS,"Writing pmatrix to " + output_str);
			CSconesIO::writeCMatrix(output_str, scones.getCMatrix(),scones.getSettings());
			logging(GIN_WARNING,"Finished in " + StringHelper::to_string<clock_t>((clock()-begin)/CLOCKS_PER_SEC) + " sec\n");

		} else {
			logging(GIN_STATUS,"Running Scones for phenotype: " + data.phenotype_names[i]);
            CScones scones(tmpData.Y.col(i),tmpData.X,tmpData.network,settings);
			scones.test_associations();

            if (debug){
                VectorXd terms = scones.getObjectiveFunctionTerms(lambda,eta);
                logging(GIN_DEBUG, "association = " + StringHelper::to_string(terms(0)));
                logging(GIN_DEBUG, "connectivity = " + StringHelper::to_string(terms(1)));
                logging(GIN_DEBUG, "sparsity = " + StringHelper::to_string(terms(2)));
                VectorXd indicator = scones.getIndicatorVector();
                VectorXd skats = scones.getScoreStatistic();
                logging(GIN_DEBUG, "Chosen SNPs\nid\tchr\tpos\tskat");
                for(uint j=0; j<indicator.rows();j++) {
                    if(indicator(j)>=1) {
                        logging(GIN_DEBUG, StringHelper::to_string(tmpData.snp_identifiers[j]) + "\t" +
                                       StringHelper::to_string(tmpData.chromosomes[j]) + "\t" +
                                       StringHelper::to_string(tmpData.positions[j]) + "\t" +
                                       StringHelper::to_string(skats(j)));
                    }
                }
            }
			logging(GIN_WARNING,"Finished in " + StringHelper::to_string<float64>(float64(clock()-begin)/CLOCKS_PER_SEC) + " sec\n");

			begin = clock();
            string output_str = outfolder_str + "/" + data.phenotype_names[i] + ".scones.pmatrix.txt";
			logging(GIN_STATUS,"Writing pmatrix to " + output_str);
			CSconesIO::writeCMatrix(output_str, scones.getCMatrix(),scones.getSettings());

			if (debug){
                VectorXd terms = scones.getObjectiveFunctionTerms(lambda,eta);
                VectorXd skat = scones.getScoreStatistic();
				output_str = outfolder_str + "/" + data.phenotype_names[i] + ".scones.out.ext.txt";
				logging(GIN_STATUS,"Writing extended output to " + output_str);
                CSconesIO::writeOutput(output_str, tmpData, scones.getIndicatorVector(),scones.getBestLambda(),scones.getBestEta(),terms,skat);

			} else {
                output_str = outfolder_str + "/" + data.phenotype_names[i] + ".scones.out.txt";
                logging(GIN_STATUS,"Writing output to " + output_str);
                CSconesIO::writeOutput(output_str, tmpData, scones.getIndicatorVector(),scones.getBestLambda(),scones.getBestEta());
            }
            logging(GIN_WARNING,"Finished in " + StringHelper::to_string<clock_t>((clock()-begin)/CLOCKS_PER_SEC) + " sec\n");
            
        }

		//logging(GIN_STATUS,"Indicator Vector:");
		//logging(GIN_INFO,scones.getIndicatorVector().transpose());
		//logging(GIN_INFO,"Best Eta: " + StringHelper::to_string<float64>(scones.getBestEta()));
		//logging(GIN_INFO,"Best Lambda: " + StringHelper::to_string<float64>(scones.getBestLambda()));
	}

	logging(GIN_STATUS,"Finished all computations in " + StringHelper::to_string<float64>((clock()-total)/CLOCKS_PER_SEC) + " sec\n");

	return 0;
}
