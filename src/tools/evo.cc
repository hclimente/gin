//
// Created by hclimente on 25/07/2017.
//

#include "gin/feature_selection/shake.h"
#include "gin/settings.h"

int main(int argc, char* argv[]) {

	Settings s(argc, argv);

	if(s.flag()) {
		s.printHelp();
		return 1;
	}

	Shake shake = Shake();

	// read the data

	shake.readGWAS(s.pedBasename(), s.encoding());
	shake.readNetwork(s.networkFilename());

	// Preprocess the data
	// shake.filterMAF(maf);
	// shake.adjustPC(pcs);

	VectorXd etas(4);
	etas << 0, 1, 2, 3;
	VectorXd lambdas(4);
	lambdas << 0, 1, 2, 3;

	shake.searchHyperparameters(10, s.modelScore(), s.associationScore());
	shake.selectSnps();

	shake.writeResults(s.output());
	return 0;

}