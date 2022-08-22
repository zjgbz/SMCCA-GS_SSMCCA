# (Supervised) Sparse Multiple Canonical Correlation Analysis 

This package provides an implementation of our sparse multiple canonical correlation analysis (SMCCA-GS) and supervised sparse multiple canonical correlation analysis (SSMCCA) methods and examples on running these two functions.

If this source code or accompanying files are helpful for your research, please cite the following publication:

## Files in this repository

All functions are in ```function_SMCCA-GS_SSMCCA.R```.


```SMCCA-GS_SSMCCA_Example.nb.html``` shows a practice of two main functions ```MultiCCA_GS``` for SMCCA-GS and ```sup_MultiCCA_GS``` for SSMCCA, and can be viewed [here](https://zjgbz.github.io/SMCCA-GS_SSMCCA/SMCCA-GS_SSMCCA_Example.nb.html).


```x1.tsv, x2.tsv, x3.tsv, y.tsv``` are data file randomly simulated data files used in the practice of ```MultiCCA_GS``` and ```sup_MultiCCA_GS```.


```.gitignore``` is the file defining which files should not be staged.

## Setup

Our functions have been tested on R 4.1.0 and R 4.1.3. All functions are self-included in ```function_SMCCA-GS_SSMCCA.R```

## Description of Main Functions

SMCCA-GS is used for finding the maximal sum correlation across all input assays without involving any outcomes, and is implemented in the function ```MultiCCA_GS```. SSMCCA is used for finding the maximal sum correlation across all input assays with considering outcomes, and is implemented in the function ```sup_MultiCCA_GS```.

### ```MultiCCA_GS```

#### Usage

````
MultiCCA_GS(moData, update_type="nores", opt_num=4, ncomponents=1, nperms=10, 
            niter_perm=3, niter=25, cca_seed=42)
````

#### Arguments

| Arguments   | Descriptions                                                                                                                                                                                                                                                                                                                                    |
|-------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| moData      | contains a list of dataframe, each of which is an assay. Each dataframe would be automatically standardized in the ```MultiCCA_GS``` by R function ```scale(T, T)```. The variance of each column of each dataframe should be non-zero, i.e., all-zero columns should be removed.                                                               |
| update_type | (optional) this is string which can either take the values "nores" or "res". "nores" means that updating weights follows the equation (17) in Witten (2009) (PMID: 19572827). "res" means that updating weights follows the appendix B in Witten (2009) (PMID: 19572827). The default value is "nores".                                         |
| opt_num     | (optional) numerical value or string "full". Numerical value means that ```MultiCCA_GS``` will re-calculate penalties for the first ```opt_num``` of CVs. If ```opt_num="full"``` or ```opt_num``` is larger than the number of CVs assigned, ```MultiCCA_GS``` will re-calculate penalties for all CVs. The default value is ```opt_num=4```.  |
| ncomponents | (optional) numerical value for the number of CVs that the user intends to calculate. The default value is ```ncomponents=1```.                                                                                                                                                                                                                  |
| nperms      | (optional) numerical value. Penalties are calculated using permutation based method. The number of permutation is determined by nperms. The default value is ```nperms=10```.                                                                                                                                                                   |
| niter_perm  | (optional) numerical value. Weights are updating in an iterative way and weights updating steps would happen in both penalty calcuation step and CV calculation step. ```niter_perm``` defining the limit of the number of iteration steps in the weight updating step of the penalty calcuation step. The default value is ```niter_perm=3```. |
| niter       | (optional) numerical value. Similar to ```niter_perm```, but ```niter``` defines the limit of the number of iteration steps in the weight updating step of the CV calcuation step. The default value is ```niter_perm=25```.                                                                                                                    |
| cca_seed    | (optional) numerical value for reproducing the results. The default value is ```cca_seed=42```.                                                                                                                                                                                                                                                 |

#### Values

| Values    | Descriptions                                                                                                                                                                                                                                |
|-----------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| canon_var | a list of dataframes contains canonical vectors (CVs) for each assay. Each CV dataframe is in the same dimension: the number of rows is the same as the number of rows of each assay, and the number of columns is the number of CVs.       |
| weight    | a list of dataframes contains weights for each assay. Each weight dataframe has the same number of columns as the number of CVs; the number of rows of each weight dataframe is the same as the number of columns of the corresponding assay. |
| penalty   | an R dataframe lists penalty of each assay of each CV. The number of columns is the same as the number of assays, and the number of rows is the same as the number of CVs.                                                                  |

### ```sup_MultiCCA_GS```

#### Usage

````
sup_MultiCCA_GS <- function(moData, y, outcome, opt_num=4, update_type="nores", assay_name_list=NULL,
                            qt_list=NULL, ncomponents=1, nperms=10, niter_perm=25, niter=3, cca_seed=42)
````

#### Arguments

| Arguments       | Descriptions                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            |
|-----------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| moData          | contains a list of dataframe, each of which is an assay. Each dataframe would be automatically standardized in the ```sup_MultiCCA_GS``` by R function ```scale(T, T)```. The variance of each column of each dataframe should be non-zero, i.e., all-zero columns should be removed.                                                                                                                                                                                                                   |
| y               | is an R dataframe or an R vector contains outcome values, or a list of vectors contains outcomes of each assay. ```sup_MultiCCA_GS``` does not consider adjustment covariates, so users can first find the residual of the outcome regressing out the adjustment covariates, and residual is the input. For each assay, the adjustment covariates might be different, so the residual might also be different. In this case, users can save them to an R list in the order of the input list of assays. |
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
