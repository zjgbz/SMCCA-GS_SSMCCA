# Sparse Multiple Canonical Correlation Analysis 

This package provides our implementations of two extensions of canonical correlation analysis (CCA): (1) sparse multiple canonical correlation analysis with Gram-Schmidt algorithm (SMCCA-GS) and (2) supervised sparse multiple canonical correlation analysis (SSMCCA). Below, we present instructions to install and run  these two functions with examples.

If this source code or accompanying files are helpful for your research, please cite the following publication:

## Files in this repository

| Filename(s)                             | Description                                                                                                                                                                                                   |
|---------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| ```function_SMCCA-GS_SSMCCA.R```      | R file that includes source codes for all functions.                                                                                                                                                                                        |
| ```SMCCA-GS_SSMCCA_Example.nb.html``` | R Markdown where two main functions ```MultiCCA_GS``` (for SMCCA-GS) and ```sup_MultiCCA_GS``` (for SSMCCA) are applied to simulated data, and can be viewed [here](https://zjgbz.github.io/SMCCA-GS_SSMCCA/SMCCA-GS_SSMCCA_Example.nb.html). |
| ```x1.tsv, x2.tsv, x3.tsv, y.tsv```   | simulated data files used in the R Markdown above.                                                                                               |

## Setup and installation

Our functions have been tested on R 4.1.0 and R 4.1.3. All functions are self-included in ```function_SMCCA-GS_SSMCCA.R```.

To install, the following line in R suffices:

````
source("function_SMCCA-GS_SSMCCA.R")
````

## Description of main functions

SMCCA-GS is used for finding canonical vectors (CVs) with maximal sum correlation across all input assays without involving any outcomes (and thus unsupervised), and is implemented in the function ```MultiCCA_GS```. SSMCCA is the supervised version where outcomes are considered, and is implemented in the function ```sup_MultiCCA_GS```.

### ```MultiCCA_GS```

#### Usage

````
MultiCCA_GS(moData, update_type="nores", opt_num=4, ncomponents=1, nperms=10, niter=25, cca_seed=42)
````

#### Arguments

| Arguments   | Description                                                                                                                                                                                                                                                                                                                                    |
|-------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| moData      | contains a list of dataframe, each of which is an assay. Each dataframe would be automatically standardized in the ```MultiCCA_GS``` by R function ```scale(T, T)```. The variance of each column of each dataframe should be non-zero: please remove non-varying columns, e.g., all-zero columns.                                                               |
| update_type | (optional) takes values "nores" or "res". With "nores", CV weights will be updated according to equation (17) in Witten (2009) (PMID: 19572827). With "res", CV weights will be updated following appendix B in Witten (2009) (PMID: 19572827). The default value is "nores".                                         |
| opt_num     | (optional) takes numerical value or a string "full". When a numerical value is specified, ```MultiCCA_GS``` will re-calculate penalties for the first ```opt_num``` of CVs. If ```opt_num="full"``` or ```opt_num``` is larger than the number of CVs assigned, ```MultiCCA_GS``` will re-calculate penalties for all CVs. The default value is ```opt_num=4```.  |
| ncomponents | (optional) numerical value for the number of CVs that the user intends to calculate. The default value is ```ncomponents=1```.                                                                                                                                                                                                                  |
| nperms      | (optional) numerical value for the number of permutations used to calibrate penalty parameters, which are used to achieve sparsity. The default value is ```nperms=10```.                                                                                                                                                                   |
| niter       | (optional) numerical value for the number of iterations used to update CV weights. Note that CV weights updating would stop early if sum correlation across all CV pairs converges before reach the ```niter```th iteration. The default value is ```niter=25```.                                                                                                                    |
| cca_seed    | (optional) numerical value for random seed used. The default value is ```cca_seed=42```.                                                                                                                                                                                                                                                 |

#### Objects

| Object(s)    | Description                                                                                                                                                                                                                                |
|-----------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| canon_var | a list of matrices that contains canonical vectors (CVs) for each assay. The dimension of a CV matrix is $n\times p$, where $n$ is the sample size of the corresponding assay, and $p$ is the number of CVs (specified by ```ncomponents```).       |
| weight    | a list of matrices that contains weights for each assay. The dimension of a weight matrix is $M\times p$, where $p$ is the number of CVs (specified by ```ncomponents```) and $M$ is the number of features of the corresponding assay. Features here are the columns of assays. |
| penalty   | a matrix that records penalty of each assay of each CV. The dimsion of a penalty matrix is $p\times K$, where $p$ is the number of CVs (specified by ```ncomponents```) and $K$ is the number of assays.                                                                  |

### ```sup_MultiCCA_GS```

#### Usage

````
sup_MultiCCA_GS(moData, y, outcome, update_type="nores", opt_num=4, assay_name_list=NULL,
                qt_list=NULL, ncomponents=1, nperms=10, niter=25, cca_seed=42)
````

#### Arguments

| Arguments       | Descriptions                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            |
|-----------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| moData          | contains a list of dataframe, each of which is an assay. Each dataframe would be automatically standardized in the ```sup_MultiCCA_GS``` by R function ```scale(T, T)```. The variance of each column of each dataframe should be non-zero, i.e., all-zero columns should be removed.                                                                                                                                                                                                                   |
| y               | is a dataframe or a vector that contains outcome values, or a list of vectors that contains outcomes of each assay. ```sup_MultiCCA_GS``` does not consider covariate adjustment explicitly. Users can adjust covariates by first obtaining residuals where covariates are regressed out of the original outcome, and then feeding residuals as y in this ```sup_MultiCCA_GS```. For each assay, the adjustment covariates might be different, so the residual might also be different. In this case, users can save them to an R list in the order of the input list of assays. |
| outcome         | is a string indicating outcome type, and can take the values of "survival", "multiclass", or "quantitative".                                                                                                                                                                                                                                                                                                                                                                                            |
| update_type     | (optional) takes values "nores" or "res". With "nores", CV weights will be updated according to equation (17) in Witten (2009) (PMID: 19572827). With "res", CV weights will be updated following appendix B in Witten (2009) (PMID: 19572827). The default value is "nores".                                                                                                                                                                                                 |
| opt_num         | (optional) takes numerical value or a string "full". When a numerical value is specified, ```sup_MultiCCA_GS``` will re-calculate penalties for the first ```opt_num``` of CVs. If ```opt_num="full"``` or ```opt_num``` is larger than the number of CVs assigned, ```sup_MultiCCA_GS``` will re-calculate penalties for all CVs. The default value is ```opt_num=4```.                                                                                                                                                  |
| qt_list         | (optional) a numerical vector that contains thresholds of keeping features for each assay. The first step of ```sup_MultiCCA_GS``` is selecting features associated with the outcomes, and ```qt_list``` provides the thresholds for each assay. The default value is 0.8 for each assay which means that top 0.2 correlated features would be kept.                                                                                                                                                                                                                    |
| ncomponents     | (optional) numerical value for the number of CVs that the user intends to calculate. The default value is ```ncomponents=1```.                                                                                                                                                                                                                                                                                                                                                                          |
| nperms          | (optional) numerical value for the number of permutations used to calibrate penalty parameters, which are used to achieve sparsity. The default value is ```nperms=10```.                                                                                                                                                                                                                                                                                                                           |
| niter           | (optional) numerical value for the number of iterations used to update CV weights. Note that CV weights updating would stop early if sum correlation across all CV pairs converges before reach the ```niter```th iteration. The default value is ```niter=25```.                                                                                                                                                                                                                                                                            |
| cca_seed        | (optional) numerical value for random seed used. The default value is ```cca_seed=42```.                                                                                                                                                                                                                                                                                                                                                                                                         |

#### Objects

| Object(s)          | Description                                                                                                                                                                                                                                                                                  |
|-----------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| canon_var       | a list of matrices that contains canonical vectors (CVs) for each assay. The dimension of a CV matrix is $n\times p$, where $n$ is the sample size of the corresponding assay, and $p$ is the number of CVs (specified by ```ncomponents```).                                                         |
| weight          | a list of matrices that contains weights for each assay. The dimension of a weight matrix is $M\times p$, where $p$ is the number of CVs (specified by ```ncomponents```) and $M$ is the number of features of the corresponding assay.                                                 |
| weight_select   | a list of matrices that contains weights for features kept for each assay. Each ```weight_selelct``` matrix has dimension $M'\times p$ where $M'$ is the number of kept features of the corresponding assay, and $p$ is the number of CVs (specified by ```ncomponents```). |
| feature_dropped | a list of vectors that contains features dropped in the feature selection step. Each vector is for one input assay.                                                                                                                                                                                                                    |
| penalty         | a matrix that records penalty of each assay of each CV. The dimsion of a penalty matrix is $p\times K$, where $p$ is the number of CVs (specified by ```ncomponents```) and $K$ is the number of assays.                                                                                                                    |

## License

This project is licensed under the GNU General Public License v3.0 - see the LICENSE.md file for details.
