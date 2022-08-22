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
| weight    | a list of dataframes contains weights for each assay. Each weight dataframe has the same number of columns: the number of CVs; the number of rows of each weight dataframe is the same as the number of columns of the corresponding assay. |
| penalty   | an R dataframe lists penalty of each assay of each CV. The number of columns is the same as the number of assays, and the number of rows is the same as the number of CVs.                                                                  |

### ```sup_MultiCCA_GS```

