# gin



[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.824641.svg)](https://doi.org/10.5281/zenodo.824641)
[![Build Status](https://travis-ci.org/hclimente/gin.svg?branch=master)](https://travis-ci.org/hclimente/gin)

gin (GWAS Incorporating Networks) is a software framework aimed at improving biomarker discovery on genotyping data using a priori information, namely networks. It is the successor of [SConES](http://bioinformatics.oxfordjournals.org/content/29/13/i171.short), the network guided multi-locus mapping method. It includes two executables (the original `scones` and `shake`, its extended version) as well as the `gin` library, ready to be used by other software, like [`martini`](https://github.com/hclimente/martini/).

## Installation

gin requires [CMake](https://cmake.org/download/) >= 3.2 to compile. To install, simply do

```
git clone --recursive https://github.com/hclimente/gin.git
gin/install_gin.sh
```

This will install `gin`, `scones` and `shake` in gin/build. If you prefer another installation path, add it is as first argument eg `gin/install_gin.sh /usr/local`.
 
## Usage

You can analyze your GWAS data from the command line executables. If you wish to use R for your analysis, please refer to the R interface [`martini`](https://github.com/hclimente/martini/). The files used in these examples are available in `test/data/case1`.

### Shake

This command is equivalent to running SConES:

```
shake --ped genotype --pheno phenotype.txt --net network.txt --depth 1
```

(See all available commands with `shake --help`.)


### SConES

This is an example of how to run SConES:

```
scones genotype phenotype.txt network.txt 0.05 . additive 0
```

The arguments are (in order):

- The prefix of the PED/MAP files, the phenotype and the network files.
- The minor allele frequency to filter.
- The output directory.
- The genetic model.
- The number of main principal components to be removed.

## Credits

gin is based on [easyGWAS](http://easygwas.ethz.ch), a C/C++ framework for computing genome-wide association studies and meta-analysis developed by [dominikgrimm](https://github.com/dominikgrimm). easyGWAS includes several standard methods for performing GWAS, such as linear regression, logistic regression and popular linear mixed models ([EMMAX](http://www.nature.com/ng/journal/v42/n4/abs/ng.548.html), [FaSTLMM](http://www.nature.com/nmeth/journal/v8/n10/abs/nmeth.1681.html)) to also account for population stratification.
