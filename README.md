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
MultiCCA_GS(moData, update_type="nores", opt_num=4, ncomponents=1, nperms=10, 
            niter_perm=3, niter=25, cca_seed=42)
````

#### Arguments

| Arguments   | Description                                                                                                                                                                                                                                                                                                                                    |
|-------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| moData      | contains a list of dataframe, each of which is an assay. Each dataframe would be automatically standardized in the ```MultiCCA_GS``` by R function ```scale(T, T)```. The variance of each column of each dataframe should be non-zero: please remove non-varying columns, e.g., all-zero columns.                                                               |
| update_type | (optional) takes values "nores" or "res". With "nores", CV weights will be updated according to equation (17) in Witten (2009) (PMID: 19572827). With "res", CV weights will be updated following appendix B in Witten (2009) (PMID: 19572827). The default value is "nores".                                         |
| opt_num     | (optional) takes numerical value or a string "full". When a numerical value is specified, ```MultiCCA_GS``` will re-calculate penalties for the first ```opt_num``` of CVs. If ```opt_num="full"``` or ```opt_num``` is larger than the number of CVs assigned, ```MultiCCA_GS``` will re-calculate penalties for all CVs. The default value is ```opt_num=4```.  |
| ncomponents | (optional) numerical value for the number of CVs that the user intends to calculate. The default value is ```ncomponents=1```.                                                                                                                                                                                                                  |
| nperms      | (optional) numerical value for the number of permutations used to calibrate penalty parameters, which are used to achieve sparsity. The default value is ```nperms=10```.                                                                                                                                                                   |
| niter_perm  | (optional) numerical value. Weights are updating in an iterative way and weights updating steps would happen in both penalty calcuation step and CV calculation step. ```niter_perm``` defining the limit of the number of iteration steps in the weight updating step of the penalty calcuation step. The default value is ```niter_perm=3```. |
| niter       | (optional) numerical value for the number of iterations used to update CV weights in the CV calculation step. The default value is ```niter_perm=25```.                                                                                                                    |
| cca_seed    | (optional) numerical value for random seed used. The default value is ```cca_seed=42```.                                                                                                                                                                                                                                                 |

#### Objects

| Object(s)    | Description                                                                                                                                                                                                                                |
|-----------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| canon_var | a list of R matrices that contains canonical vectors (CVs) for each assay. The dimension of a CV matrix is $n\times p$, where $n$ is the sample size of the corresponding assay, and $p$ is the number of CVs (specified by ```ncomponents```).       |
| weight    | a list of R matrices that contains weights for each assay. The dimension of a weight matrix is $p\times M$, where $p$ is the number of CVs (specified by ```ncomponents```) and $M$ is the number of features of the corresponding assay. Features here are the columns of assays. |
| penalty   | an R matrix that records penalty of each assay of each CV. The dimsion of a penalty matrix is $p\times K$, where $p$ is the number of CVs (specified by ```ncomponents```) and $K$ is the number of assays.                                                                  |

### ```sup_MultiCCA_GS```

#### Usage

````
sup_MultiCCA_GS(moData, y, outcome, opt_num=4, update_type="nores", assay_name_list=NULL,
                qt_list=NULL, ncomponents=1, nperms=10, niter_perm=25, niter=3, cca_seed=42)
````

#### Arguments

| Arguments       | Descriptions                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            |
|-----------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| moData          | contains a list of dataframe, each of which is an assay. Each dataframe would be automatically standardized in the ```sup_MultiCCA_GS``` by R function ```scale(T, T)```. The variance of each column of each dataframe should be non-zero, i.e., all-zero columns should be removed.                                                                                                                                                                                                                   |
| y               | is a dataframe or a vector that contains outcome values, or a list of vectors that contains outcomes of each assay. ```sup_MultiCCA_GS``` does not consider covariate adjustment explicitly. Users can adjust covariates by first obtaining residuals where covariates are regressed out of the original outcome, and then feeding residuals as y in this ```sup_MultiCCA_GS```. For each assay, the adjustment covariates might be different, so the residual might also be different. In this case, users can save them to an R list in the order of the input list of assays. |
| outcome         | is a string indicating outcome type, and can take the values of "survival", "multiclass", or "quantitative".                                                                                                                                                                                                                                                                                                                                                                                            |
| opt_num         | (optional) numerical value or string "full". Numerical value means that ```sup_MultiCCA_GS``` will re-calculate penalties for the first ```opt_num``` of CVs. If ```opt_num="full"``` or ```opt_num``` is larger than the number of CVs assigned, ```sup_MultiCCA_GS``` will re-calculate penalties for all CVs. The default value is ```opt_num=4```.                                                                                                                                                  |
| update_type     | (optional) this is string which can either take the values "nores" or "res". "nores" means that updating weights follows the equation (17) in Witten (2009) (PMID: 19572827). "res" means that updating weights follows the appendix B in Witten (2009) (PMID: 19572827). The default value is "nores".                                                                                                                                                                                                 |
| assay_name_list | (optional) this is a vector of strings contains the names of each assay. The default value is ```assay_name_list=NULL```. If users assign vector of assay names, the list of features kept in the analysis would have corresponding names. If ```assay_name_list=NULL```, the names of features dropped in the analysis would be 1, 2, 3, ..., which is the same as the index of the list of input assays. Features here are the columns of assays.                                                        |
| qt_list         | (optional) a numerical vector contains thresholds of keeping features for each assay. The first step of ```sup_MultiCCA_GS``` is selecting features associated with the outcomes, and ```qt_list``` provides the thresholds for each assay. The default value is 0.8 for each assay.                                                                                                                                                                                                                    |
| ncomponents     | (optional) numerical value for the number of CVs that the user intends to calculate. The default value is ```ncomponents=1```.                                                                                                                                                                                                                                                                                                                                                                          |
| nperms          | (optional) numerical value. Penalties are calculated using permutation based method. The number of permutation is determined by nperms. The default value is ```nperms=10```.                                                                                                                                                                                                                                                                                                                           |
| niter_perm      | (optional) numerical value. Weights are updating in an iterative way and weights updating steps would happen in both penalty calcuation step and CV calculation step. ```niter_perm``` defining the limit of the number of iteration steps in the weight updating step of the penalty calcuation step. The default value is ```niter_perm=3```.                                                                                                                                                         |
| niter           | (optional) numerical value. Similar to ```niter_perm```, but ```niter``` defines the limit of the number of iteration steps in the weight updating step of the CV calcuation step. The default value is ```niter_perm=25```.                                                                                                                                                                                                                                                                            |
| cca_seed        | (optional) numerical value for reproducing the results. The default value is ```cca_seed=42```.                                                                                                                                                                                                                                                                                                                                                                                                         |

#### Values

| Values          | Descriptions                                                                                                                                                                                                                                                                                  |
|-----------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| canon_var       | a list of dataframes contains canonical vectors (CVs) for each assay. Each CV dataframe is in the same dimension: the number of rows is the same as the number of rows of each assay, and the number of columns is the number of CVs.                                                         |
| weight          | a list of dataframes contains weights for each assay. Each weight dataframe has the same number of columns as the number of CVs; the number of rows of each weight dataframe is the same as the number of columns of the corresponding assay.                                                 |
| weight_select   | a list of dataframes contains weights of kept features for each assay. The number of columns of each output dataframe is the number of kept features of the corresponding assay; the number of rows of each output dataframe is the same as the number of columns of the corresponding assay. |
| feature_dropped | a list of vectors contains features dropped in the feature selection step.                                                                                                                                                                                                                    |
| penalty         | an R dataframe lists penalty of each assay of each CV. The number of columns is the same as the number of assays, and the number of rows is the same as the number of CVs.                                                                                                                    |

## License

This project is licensed under the GNU General Public License v3.0 - see the LICENSE.md file for details.
