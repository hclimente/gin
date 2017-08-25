//
// Created by hclimente on 25/07/2017.
//

#include "gin/feature_selection/shake.h"
#include "gin/settings.h"

int main(int argc, char* argv[]) {

	Settings s(argc, argv);

	if(s.error()) {
		s.printHelp();
		return 1;
	}

	Shake experiment = Shake();
	experiment.setDebug(s.debug());

	// read the data

	experiment.readGWAS(s.pedBasename(), s.encoding());
	experiment.readNetwork(s.networkFilename());

	// Preprocess the data
	// experiment.filterMAF(maf);
	// experiment.adjustPC(pcs);

	VectorXd etas(4);
	etas << 0, 1, 2, 3;
	VectorXd lambdas(4);
	lambdas << 0, 1, 2, 3;

	if (s.etas().rows() == 0) {
		experiment.selectHyperparameters(10, s.modelScore(), s.associationScore());
	} else {
		experiment.selectHyperparameters(10, s.modelScore(), s.associationScore(), s.etas(), s.lambdas());
	}

	experiment.selectSNPs();

	experiment.writeResults(s.output());
	return 0;

}