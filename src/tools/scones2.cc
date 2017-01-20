#include <iostream>
#include "boost/program_options.hpp"

#include "CEasyGWAS/gwas/CScones.h"
#include "CEasyGWAS/io/CPlinkParser.h"
#include "CEasyGWAS/kernel/CKernels.h"
#include "CEasyGWAS/io/CSconesIO.h"
#include "CEasyGWAS/stats/CStats.h"
#include "CEasyGWAS/utils/CMatrixHelper.h"

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
		("maf,m", po::value<double>()->default_value(0.05), "Float minor allele frequency filter.")
		("outdir,o", po::value<string>()->default_value("."), "Output directory.")
		("encoding,s", po::value<string>()->default_value("additive"), "snp_encoding.")
		("pc,c", po::value<int>()->default_value(0), "PC.")
		("lambda,l", po::value<double>()->default_value(-1), "Lambda parameter.")
		("eta,e", po::value<double>()->default_value(-1), "Eta parameter.")
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
	double lambda = vm["lambda"].as<double>();
	double eta = vm["eta"].as<double>();
    bool debug = vm["debug"].as<bool>();

	uint encoding = 0;
	if(snp_encoding=="additive") encoding = 0;
	else if(snp_encoding=="recessive") encoding = 1;
	else if(snp_encoding=="dominant") encoding = 2;
	else if(snp_encoding=="overdominant") encoding = 3;
	else {
		logging(ERROR,"Encoding does not exist!");
		exit(-1);
	}

    GWASData data;

	float64 begin;
	begin = clock();
	logging(STATUS,"Reading Genotype file...");
	CPlinkParser::readPEDFile(genotype_str + ".ped", &data);
	logging(INFO,"Number of SNPs: " + StringHelper::to_string<uint64>(data.n_snps));
	logging(INFO,"Number of Samples: " + StringHelper::to_string<uint64>(data.n_samples));
	logging(WARNING,"Finished in " + StringHelper::to_string<float64>(float64(clock()-begin)/CLOCKS_PER_SEC) + " sec\n");

	begin = clock();
	logging(STATUS,"Reading Mapping file...");
	CPlinkParser::readMAPFile(genotype_str + ".map", &data);
	logging(WARNING,"Finished in " + StringHelper::to_string<float64>(float64(clock()-begin)/CLOCKS_PER_SEC) + " sec\n");


	begin = clock();
	logging(STATUS,"Reading Phenotype file...");
	CPlinkParser::readPhenotypeFile(phenotype_str,&data);
	logging(WARNING,"Finished in " + StringHelper::to_string<float64>(float64(clock()-begin)/CLOCKS_PER_SEC) + " sec\n");


	begin = clock();
	logging(STATUS,"Encoding SNP data...");
	CGWASDataHelper::encodeHeterozygousData(&data,encoding);
	logging(WARNING,"Finished in " + StringHelper::to_string<float64>(float64(clock()-begin)/CLOCKS_PER_SEC) + " sec\n");

	/*
	begin = clock();
	logging(STATUS,"Filter unique SNPs...");
	CGWASDataHelper::filterUniqueSNPs(&data);
	logging(WARNING,"Finished in " + StringHelper::to_string<float64>((clock()-begin)/CLOCKS_PER_SEC) + " sec\n");
	*/

	begin = clock();
	logging(STATUS,"Filter SNPs by MAF...");
	CGWASDataHelper::filterSNPsByMAF(&data,maf);
	logging(INFO,"Number of SNPs (after MAF): " + StringHelper::to_string<uint64>(data.n_snps));
	logging(WARNING,"Finished in " + StringHelper::to_string<float64>(float64(clock()-begin)/CLOCKS_PER_SEC) + " sec\n");

	begin = clock();
	logging(STATUS,"Loading and filtering network file...");
	CSconesIO::readSparseNetworkFile(network_str,&data);
	logging(WARNING,"Finished in " + StringHelper::to_string<float64>(float64(clock()-begin)/CLOCKS_PER_SEC) + " sec\n");


	for(uint i=0; i<data.phenotype_names.size(); i++) {

		begin = clock();
		logging(STATUS,"Remove samples with missing values...");
		GWASData tmpData = CGWASDataHelper::removeSamples4MissingData(data,i);
		if(tmpData.n_samples!=data.n_samples) {
			logging(INFO,"#Samples removed: " + StringHelper::to_string<uint64>(data.n_samples-tmpData.n_samples));
		}
		logging(WARNING,"Finished in " + StringHelper::to_string<float64>(float64(clock()-begin)/CLOCKS_PER_SEC) + " sec\n");
		logging(WARNING,"Total Data: " + StringHelper::to_string<float64>(float64(clock()-total_data)/CLOCKS_PER_SEC) + " sec\n");

		//OPTIONAL SET DIFFERENT SConES Settings
		begin = clock();
		//COMPUTE PRINCIPLE COMPONENTS
		if (pcs>0) {
			begin = clock();
			logging(STATUS,"Computing Realized Relationship Kernel and Principle Components...");
			data.K = CKernels::realizedRelationshipKernel(data.X);
			MatrixXd PCs = CStats::principle_components(data.K);
			PCs = sliceColsMatrix(PCs,VectorXd::LinSpaced(pcs,0,pcs-1));
			logging(WARNING,"Finished in " + StringHelper::to_string<clock_t>((clock()-begin)/CLOCKS_PER_SEC) + " sec\n");

			logging(STATUS,"Running Scones for phenotype: " + data.phenotype_names[i]);
			CSconesSettings settings;
			CScones scones(tmpData.Y.col(i),tmpData.X,tmpData.network,PCs,settings);
			if (lambda == -1 && eta == -1) scones.test_associations();
			else scones.test_associations(lambda, eta);
			logging(WARNING,"Finished in " + StringHelper::to_string<float64>(float64(clock()-begin)/CLOCKS_PER_SEC) + " sec\n");

			begin = clock();
			string output_str = outfolder_str + "/" + data.phenotype_names[i] + ".scones.out.txt";
			logging(STATUS,"Writing output to " + output_str);
            if (lambda == -1 && eta == -1) CSconesIO::writeOutput(output_str, tmpData, scones.getIndicatorVector(),scones.getBestLambda(),scones.getBestEta());
            else CSconesIO::writeOutput(output_str, tmpData, scones.getIndicatorVector(),lambda,eta);
			output_str = outfolder_str + "/" + data.phenotype_names[i] + ".scones.pmatrix.txt";
			logging(STATUS,"Writing pmatrix to " + output_str);
			CSconesIO::writeCMatrix(output_str, scones.getCMatrix(),scones.getSettings());
			logging(WARNING,"Finished in " + StringHelper::to_string<clock_t>((clock()-begin)/CLOCKS_PER_SEC) + " sec\n");

		} else {
			logging(STATUS,"Running Scones for phenotype: " + data.phenotype_names[i]);
			CScones scones(tmpData.Y.col(i),tmpData.X,tmpData.network);
			if (lambda == -1 && eta == -1) scones.test_associations();
			else {
                scones.test_associations(lambda, eta);

                if (debug){
                    VectorXd terms = scones.getObjectiveFunctionTerms(lambda,eta);
                    logging(DEBUG, "association = " + StringHelper::to_string(terms(0)));
                    logging(DEBUG, "connectivity = " + StringHelper::to_string(terms(1)));
                    logging(DEBUG, "sparsity = " + StringHelper::to_string(terms(2)));
                    VectorXd indicator = scones.getIndicatorVector();
                    VectorXd skats = scones.getScoreStatistic();
                    logging(DEBUG, "Chosen SNPs\nid\tchr\tpos\tskat");
                    for(uint j=0; j<indicator.rows();j++) {
                        if(indicator(j)>=1) {
                            logging(DEBUG, StringHelper::to_string(tmpData.snp_identifiers[j]) + "\t" +
									StringHelper::to_string(tmpData.chromosomes[j]) + "\t" +
									StringHelper::to_string(tmpData.positions[j]) + "\t" +
									StringHelper::to_string(skats(j)));
                        }
                    }
                }
            }
			logging(WARNING,"Finished in " + StringHelper::to_string<float64>(float64(clock()-begin)/CLOCKS_PER_SEC) + " sec\n");

			begin = clock();
			string output_str = outfolder_str + "/" + data.phenotype_names[i] + ".scones.out.txt";
			logging(STATUS,"Writing output to " + output_str);
			if (lambda == -1 && eta == -1) CSconesIO::writeOutput(output_str, tmpData, scones.getIndicatorVector(),scones.getBestLambda(),scones.getBestEta());
			else CSconesIO::writeOutput(output_str, tmpData, scones.getIndicatorVector(),lambda,eta);
			output_str = outfolder_str + "/" + data.phenotype_names[i] + ".scones.pmatrix.txt";
			logging(STATUS,"Writing pmatrix to " + output_str);
			CSconesIO::writeCMatrix(output_str, scones.getCMatrix(),scones.getSettings());
			logging(WARNING,"Finished in " + StringHelper::to_string<clock_t>((clock()-begin)/CLOCKS_PER_SEC) + " sec\n");

			if (debug){
                VectorXd terms = scones.getObjectiveFunctionTerms(lambda,eta);
                VectorXd skat = scones.getScoreStatistic();
				output_str = outfolder_str + "/" + data.phenotype_names[i] + ".scones.out.ext.txt";
				logging(STATUS,"Writing extended output to " + output_str);
                if (lambda == -1 && eta == -1) CSconesIO::writeOutput(output_str, tmpData, scones.getIndicatorVector(),scones.getBestLambda(),scones.getBestEta(),terms,skat);
                else CSconesIO::writeOutput(output_str, tmpData, scones.getIndicatorVector(),lambda,eta,terms,skat);
				logging(WARNING,"Finished in " + StringHelper::to_string<clock_t>((clock()-begin)/CLOCKS_PER_SEC) + " sec\n");

                output_str = outfolder_str + "/" + data.phenotype_names[i] + ".scones.L.txt";
                logging(STATUS,"Writing Laplacian matrix to " + output_str);
                CSconesIO::writeLaplacianMatrix(output_str, tmpData, tmpData.network);
                logging(WARNING,"Finished in " + StringHelper::to_string<clock_t>((clock()-begin)/CLOCKS_PER_SEC) + " sec\n");

			}


		}

		//logging(STATUS,"Indicator Vector:");
		//logging(INFO,scones.getIndicatorVector().transpose());
		//logging(INFO,"Best Eta: " + StringHelper::to_string<float64>(scones.getBestEta()));
		//logging(INFO,"Best Lambda: " + StringHelper::to_string<float64>(scones.getBestLambda()));
	}

	logging(STATUS,"Finished all computations in " + StringHelper::to_string<float64>((clock()-total)/CLOCKS_PER_SEC) + " sec\n");

	return 0;
}
