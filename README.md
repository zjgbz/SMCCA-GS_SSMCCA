# (Supervised) Sparse Multiple Canonical Correlation Analysis 

This package provides an implementation of our sparse multiple canonical correlation analysis (SMCCA-GS) and supervised sparse multiple canonical correlation analysis (SSMCCA) methods and examples on running these two functions.

A practice of these two function is [here](https://zjgbz.github.io/SMCCA-GS_SSMCCA/SMCCA-GS_SSMCCA_Example.nb.html).

If this source code or accompanying files are helpful for your research, please cite the following publication:

## Setup

Our functions have been tested on R 4.1.0 and R 4.1.3. All functions are self-included in ```function_SMCCA-GS_SSMCCA.R```

## Description of Main Functions

SMCCA-GS is used for finding the maximal sum correlation across all input assays without involving any outcomes, and is implemented in the function ```MultiCCA_GS```. SSMCCA is used for finding the maximal sum correlation across all input assays with considering outcomes, and is implemented in the function ```sup_MultiCCA_GS```.

### ```MultiCCA_GS```
SMCCA-GS is used for finding the maximal sum correlation across all input assays without involving any outcomes, and is implemented in the function ```MultiCCA_GS```.


## Abbreviations:
PMA
SMCCA
SSCCA
GS
SMCCA-GS
SSMCCA

Multi-passes version of supervised spare CCA based on the modified PMA package.
