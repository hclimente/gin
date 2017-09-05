#include "gin/stats/univariate_association.h"

class MockUnivariateAssociation : public UnivariateAssociation {

public:
	MOCK_METHOD0(computeSKAT, VectorXd());
	MOCK_METHOD1(computeSKAT, VectorXd(VectorXd));
	MOCK_METHOD0(computeChi2, VectorXd());
	MOCK_METHOD1(computeTrendTest, VectorXd(std::string const&));

};
