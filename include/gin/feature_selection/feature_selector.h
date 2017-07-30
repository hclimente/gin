//
// Created by hclimente on 23/07/2017.
//

#ifndef GIN_FEATURE_SELECTION_H
#define GIN_FEATURE_SELECTION_H

#include "gin/globals.h"

class FeatureSelector {
public:
	FeatureSelector(long);

	VectorXd selected();

	void setSelected(VectorXd);

protected:

	VectorXd __selected;
	long __n_features;

};

#endif //GIN_FEATURE_SELECTION_H
