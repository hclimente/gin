#include "gtest/gtest.h"
#include "gin/gwas/CScones.h"
#include "gin/io/CSconesIO.h"
#include "gin/io/CPlinkParser.h"
#include "gin/globals.h"

struct CSconesInitialSettings{

    float64 eta;
    float64 expected_eta;
    float64 lambda;
    float64 expected_lambda;
    int griddepth;
    uint test_statistic;
    float64 expected_association;
    float64 expected_connectivity;
    float64 expected_sparsity;
    string path_prefix;
    int selection_criterion;
    vector<int> expected_causal_SNPs;

};

#ifdef __linux__

CSconesInitialSettings skat_grid_consistency = CSconesInitialSettings {
        -1, // eta
        16681.00537, // expected_eta
        -1, // lambda
        278.25594, // expected_lambda
        1, // griddepth
        SKAT, // test_statistic
        724379, // expected_association
        44799.20, // expected_connectivity
        166810, // expected_sparsity
        "test/data/case1/", // path_prefix
        CONSISTENCY, // selection_criterion
        {676, 679, 680, 682, 684, 685, 686, 690, 695, 696} // expected_causal_SNPs
};

CSconesInitialSettings skat_grid_information = CSconesInitialSettings {
        -1, // eta
        16681.00537, // expected_eta
        -1, // lambda
        2154.43469, // expected_lambda
        1, // griddepth
        SKAT, // test_statistic
        820863, // expected_association
        4308, // expected_connectivity
        433706, // expected_sparsity
        "test/data/case1/", // path_prefix
        AICc, // selection_criterion
        {676, 677, 678, 679, 680, 681, 682, 683, 684, 685, 686, 687, 688, 689, 690, 691, 692, 693, 694, 695, 696, 697, 698, 699, 700, 701} // expected_causal_SNPs
};

CSconesInitialSettings skat_fixed = CSconesInitialSettings {
        17000, // eta
        17000, // expected_eta
        300, // lambda
        300, // expected_lambda
        1, // griddepth
        SKAT, // test_statistic
        724379, // expected_association
        48300, // expected_connectivity
        170000, // expected_sparsity
        "test/data/case1/", // path_prefix
        CONSISTENCY, // selection_criterion
        {676, 679, 680, 682, 684, 685, 686, 690, 695, 696} // expected_causal_SNPs
};

CSconesInitialSettings chisq_fixed = CSconesInitialSettings {
        13.0694, // eta
        13.0694, // expected_eta
        0.218, // lambda
        0.218, // expected_lambda
        1, // griddepth
        CHISQ, // test_statistic
        517.83, // expected_association
        35.1, // expected_connectivity
        130.694, // expected_sparsity
        "test/data/case1/", // path_prefix
        CONSISTENCY, // selection_criterion
        {676, 679, 680, 682, 684, 685, 686, 690, 695, 696} // expected_causal_SNPs

};

CSconesInitialSettings chisq_grid5_information = CSconesInitialSettings {
        -1, // eta
        33435.468, // expected_eta
        -1, // lambda
        3.385, // expected_lambda
        5, // griddepth
        SKAT, // test_statistic
        625692, // expected_association
        453.695, // expected_connectivity
        234048, // expected_sparsity
        "test/data/case1/", // path_prefix
        CONSISTENCY, // selection_criterion
        {676, 680, 682, 685, 686, 690, 696} // expected_causal_SNPs

};

#endif

#ifdef __APPLE__

CSconesInitialSettings skat_fixed = CSconesInitialSettings {
        17000, // eta
        17000, // expected_eta
        300, // lambda
        300, // expected_lambda
        1, // griddepth
        SKAT, // test_statistic
        693163, // expected_association
        46200, // expected_connectivity
        153000, // expected_sparsity
        "test/data/case1/", // path_prefix
        CONSISTENCY, // selection_criterion
        {676, 679, 680, 682, 685, 686, 690, 695, 696} // expected_causal_SNPs
};

CSconesInitialSettings skat_grid_consistency = CSconesInitialSettings {
        -1, // eta
        16681.00537, // expected_eta
        -1, // lambda
        278.25594, // expected_lambda
        1, // griddepth
        SKAT, // test_statistic
        724379, // expected_association
        44799, // expected_connectivity
        166810, // expected_sparsity
        "test/data/case1/", // path_prefix
        CONSISTENCY, // selection_criterion
        {676, 679, 680, 682, 684, 685, 686, 690, 695, 696} // expected_causal_SNPs
};

CSconesInitialSettings skat_grid_information = CSconesInitialSettings {
        -1, // eta
        16681.00537, // expected_eta
        -1, // lambda
        2154.43469, // expected_lambda
        1, // griddepth
        SKAT, // test_statistic
        820863, // expected_association
        4308, // expected_connectivity
        433706, // expected_sparsity
        "test/data/case1/", // path_prefix
        AICc, // selection_criterion
        {676, 677, 678, 679, 680, 681, 682, 683, 684, 685, 686, 687, 688, 689, 690, 691, 692, 693, 694, 695, 696, 697, 698, 699, 700, 701} // expected_causal_SNPs
};

CSconesInitialSettings chisq_fixed = CSconesInitialSettings {
        13.0694, // eta
        13.0694, // expected_eta
        0.218, // lambda
        0.218, // expected_lambda
        1, // griddepth
        CHISQ, // test_statistic
        517.83, // expected_association
        35.1, // expected_connectivity
        130.694, // expected_sparsity
        "test/data/case1/", // path_prefix
        CONSISTENCY, // selection_criterion
        {676, 679, 680, 682, 684, 685, 686, 690, 695, 696} // expected_causal_SNPs

};

CSconesInitialSettings chisq_grid5_information = CSconesInitialSettings {
        -1, // eta
        33435.468, // expected_eta
        -1, // lambda
        3.385, // expected_lambda
        5, // griddepth
        SKAT, // test_statistic
        625692, // expected_association
        453.695, // expected_connectivity
        234048, // expected_sparsity
        "test/data/testing/scones/case1/", // path_prefix
        CONSISTENCY, // selection_criterion
        {676, 677, 678, 679, 680, 681, 682, 683, 684, 685, 686, 687, 688, 689, 690, 691, 692, 693, 694, 695, 696, 697, 698, 699, 700, 701} // expected_causal_SNPs

};

#endif