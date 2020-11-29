rm(list=ls())
library(PMA)
#source("/Users/jiangminzhi/Dropbox/Yun_Li_project/mesa_cca/mcca_script/CCA.R")
source("sup_MultiCCA.R")
# source("CCA_custom.R")
library(feather)

X2 = read_feather("common_exam5-pheno01-auto_log-prot_methy_BIS_log-rnaseq_RSEMv1.3.0.rsem_genes_tpm_PBMC_remove_outliers_feature_expressed.feather")
X1 = read_feather("common_exam5-pheno01-auto_log-rnaseq_RSEMv1.3.0.rsem_genes_tpm_PBMC_log-prot_methy_BIS_remove_outliers_feature_expressed.feather")
pheno = read.table("common_exam5-pheno01_pheno01_prot_methy_rnaseq_remove_outliers.tsv", sep = "\t", header = 1)
y = pheno["wbc_ncnc_bld_1"]

print(dim(X1))
X1 = X1[, colSums(X1 != 0) > 0]
print(dim(X1))
X1_norm = scale(X1, T, T)
X2_norm = scale(X2, T, T)

moData = list(X1, X2)
outcome = "quantitative"
num_components = 15

fullMCCA_1 = sup_MultiCCA.permute(moData, y, outcome)
fullMCCA2_1 = sup_MultiCCA(moData, y, outcome, penalty = fullMCCA_1$bestpenalties, ncomponents = num_components)
X1_sel = fullMCCA2_1$xlist[[1]]
X2_sel = fullMCCA2_1$xlist[[2]]
X1_sel_norm = scale(X1_sel, T, T)
X2_sel_norm = scale(X2_sel, T, T)

w1_sup_mcca = fullMCCA2_1$ws[[1]]
w2_sup_mcca = fullMCCA2_1$ws[[2]]

w1_sup_mcca_nonzero = w1_sup_mcca[rowSums(w1_sup_mcca != 0) > 0, ]
w2_sup_mcca_nonzero = w2_sup_mcca[rowSums(w2_sup_mcca != 0) > 0, ]

var1_sup_mcca = as.data.frame(as.matrix(X1_sel_norm) %*% w1_sup_mcca)
var2_sup_mcca = as.data.frame(as.matrix(X2_sel_norm) %*% w2_sup_mcca)
cor_var1_var2_sup_mcca = cor(var1_sup_mcca[, 1], var2_sup_mcca[, 1])
cor_y_var1_sup_mcca = cor(var1_sup_mcca[, 1], y)
cor_y_var2_sup_mcca = cor(var2_sup_mcca[, 1], y)

perm.out = CCA.permute(x=X1, z=X2, typex="standard",typez="standard",
                       outcome=outcome, y=y)
out <- CCA(x=X1,z=X2,typex="standard",typez="standard",K=num_components,
           penaltyx=perm.out$bestpenaltyx,penaltyz=perm.out$bestpenaltyz,
           v=perm.out$v.init, outcome=outcome, y=y)

w1_sup_cca = out$u
w2_sup_cca = out$v

w1_sup_cca_nonzero = w1_sup_cca[rowSums(w1_sup_cca != 0) > 0, ]
w2_sup_cca_nonzero = w2_sup_cca[rowSums(w2_sup_cca != 0) > 0, ]

var1_sup_cca = as.data.frame(as.matrix(X1_norm) %*% w1_sup_cca)
var2_sup_cca = as.data.frame(as.matrix(X2_norm) %*% w2_sup_cca)
cor_var1_var2_sup_cca = cor(var1_sup_cca[, 1], var2_sup_cca[, 1])
cor_y_var1_sup_cca = cor(var1_sup_cca[, 1], y)
cor_y_var2_sup_cca = cor(var2_sup_cca[, 1], y)

cor_var1_sup_cca_mcca = cor(var1_sup_cca, var1_sup_mcca)
cor_var2_sup_cca_mcca = cor(var2_sup_cca, var2_sup_mcca)

