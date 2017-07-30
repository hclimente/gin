//
// Created by hclimente on 23/07/2017.
//

//
// Created by hclimente on 21/07/2017.
//

#include "gin/feature_selection/scones.h"
#include "maxflow/maxflow.h"

FeatureSelector::FeatureSelector(long n_features) {
	__n_features = n_features;
	__selected = VectorXd::Ones(__n_features);
}

VectorXd FeatureSelector::selected() {
	return __selected;
}

void FeatureSelector::setSelected(VectorXd selected) {
	__selected = selected;
}