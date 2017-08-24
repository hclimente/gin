//
// Created by hclimente on 24/08/2017.
//

class MockGridCV : public GridCV {
public:
	MOCK_METHOD1(runFolds, void(uint));
	MOCK_METHOD2(runFolds, void(uint, CCrossValidation));
	MOCK_METHOD1(scoreModels, void(uint));
	MOCK_METHOD0(bestEta, double());
	MOCK_METHOD0(bestLambda, double());
	MOCK_METHOD0(etas, VectorXd());
	MOCK_METHOD0(lambdas, VectorXd());
	MOCK_METHOD0(grids, std::vector<Grid*>());
	MOCK_METHOD0(binary_y, bool());
	MOCK_METHOD0(aggregatedFolds, std::map<double, std::map<double, VectorXd>>());
	MOCK_METHOD2(aggregatedFolds, VectorXd(uint const& e, uint const& l));
	MOCK_METHOD0(scoredFolds, MatrixXd());
	MOCK_METHOD2(scoredFolds, double(uint const& e, uint const& l));
};