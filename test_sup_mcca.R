rm(list=ls())
#library(PMA)
#source("/Users/jiangminzhi/Dropbox/Yun_Li_project/mesa_cca/mcca_script/CCA.R")
source("sup_MultiCCA.R")
source("CCA_custom.R")

set.seed(11)
X1 = read.table("X1.tsv", sep = "\t", header = FALSE)
X2 = read.table("X2.tsv", sep = "\t", header = FALSE)
y = read.table("y.tsv", sep = "\t", header = FALSE)

X1_norm = scale(X1, T, T)
X2_norm = scale(X2, T, T)

moData = list(X1, X2)
outcome = "quantitative"
num_components = 10

# xlist = MultiCCA.Phenotype.ZeroSome(moData, y, qt = .8, cens = NULL, outcome = outcome, type = c("standard", "standard"))
# m_x = xlist$xlist_sel[[1]]
# m_z = xlist$xlist_sel[[2]]
# pheno.out <- CCAPhenotypeZeroSome(X1,X2,y,qt=.8, cens= NULL, outcome=outcome, typex="standard", typez="standard")
# x = pheno.out$x
# z = pheno.out$z
# x_effec = x[, colSums(x != 0) > 0]
# z_effec = z[, colSums(z != 0) > 0]

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
                       outcome=outcome, y=y, presel = FALSE)
out <- CCA(x=X1,z=X2,typex="standard",typez="standard",K=num_components,
           penaltyx=perm.out$bestpenaltyx,penaltyz=perm.out$bestpenaltyz,
           v=perm.out$v.init, outcome=outcome, y=y, presel = FALSE)

w1_sup_cca = out$u
w2_sup_cca = out$v

w1_sup_cca_nonzero = w1_sup_cca[rowSums(w1_sup_cca != 0) > 0, ]
w2_sup_cca_nonzero = w2_sup_cca[rowSums(w2_sup_cca != 0) > 0, ]

var1_sup_cca = as.data.frame(as.matrix(X1_norm) %*% w1_sup_cca)
var2_sup_cca = as.data.frame(as.matrix(X2_norm) %*% w2_sup_cca)
cor_var1_var2_sup_cca = cor(var1_sup_cca[, 1], var2_sup_cca[, 1])
cor_y_var1_sup_cca = cor(var1_sup_cca[, 1], y)
cor_y_var2_sup_cca = cor(var2_sup_cca[, 1], y)

# perm.out_presel = CCA.permute(x=X1_sel, z=X2_sel, typex="standard",typez="standard",
#                        outcome=outcome, y=y, presel = TRUE)
# out_presel <- CCA(x=X1_sel,z=X2_sel,typex="standard",typez="standard",K=num_components,
#            penaltyx=perm.out_presel$bestpenaltyx,penaltyz=perm.out_presel$bestpenaltyz,
#            v=perm.out_presel$v.init, outcome=outcome, y=y, presel = TRUE)

# w1_sup_cca_presel = out_presel$u
# w2_sup_cca_presel = out_presel$v

# w1_sup_cca_presel_nonzero = w1_sup_cca_presel[rowSums(w1_sup_cca_presel != 0) > 0, ]
# w2_sup_cca_presel_nonzero = w2_sup_cca_presel[rowSums(w2_sup_cca_presel != 0) > 0, ]

# var1_sup_cca_presel = as.data.frame(as.matrix(X1_sel_norm) %*% w1_sup_cca_presel)
# var2_sup_cca_presel = as.data.frame(as.matrix(X2_sel_norm) %*% w2_sup_cca_presel)
# cor_var1_var2_sup_cca_presel = cor(var1_sup_cca_presel[, 1], var2_sup_cca_presel[, 1])
# cor_y_var1_sup_cca_presel = cor(var1_sup_cca_presel[, 1], y)
# cor_y_var2_sup_cca_presel = cor(var2_sup_cca_presel[, 1], y)

# unsupervisee

# fullMCCA_2 = MultiCCA.permute(moData)
# fullMCCA2_2 = MultiCCA(moData, penalty = fullMCCA_2$bestpenalties, ncomponents = num_components)

# w1_unsup_mcca = fullMCCA2_2$ws[[1]]
# w2_unsup_mcca = fullMCCA2_2$ws[[2]]
# var1_unsup_mcca = as.data.frame(as.matrix(X1_norm) %*% w1_unsup_mcca)
# var2_unsup_mcca = as.data.frame(as.matrix(X2_norm) %*% w2_unsup_mcca)
# cor_var1_var2_unsup_mcca = cor(var1_unsup_mcca[, 1], var2_unsup_mcca[, 1])
# cor_y_var1_unsup_mcca = cor(var1_unsup_mcca[, 1], y)
# cor_y_var2_unsup_mcca = cor(var2_unsup_mcca[, 1], y)

# moData_sel = list(X1_sel, X2_sel)
# fullMCCA_3 = MultiCCA.permute(moData_sel, nperms=25)
# fullMCCA2_3 = MultiCCA(moData_sel, niter=15, penalty = fullMCCA_3$bestpenalties, ncomponents = num_components)

# w1_unsup_mcca_sel = fullMCCA2_3$ws[[1]]
# w2_unsup_mcca_sel = fullMCCA2_3$ws[[2]]
# var1_unsup_mcca_sel = as.data.frame(as.matrix(X1_sel_norm) %*% w1_unsup_mcca_sel)
# var2_unsup_mcca_sel = as.data.frame(as.matrix(X2_sel_norm) %*% w2_unsup_mcca_sel)
# cor_var1_var2_unsup_mcca_sel = cor(var1_unsup_mcca_sel[, 1], var2_unsup_mcca_sel[, 1])
# cor_y_var1_unsup_mcca_sel = cor(var1_unsup_mcca_sel[, 1], y)
# cor_y_var2_unsup_mcca_sel = cor(var2_unsup_mcca_sel[, 1], y)

cor_var1_sup_cca_mcca = cor(var1_sup_cca, var1_sup_mcca)
cor_var2_sup_cca_mcca = cor(var2_sup_cca, var2_sup_mcca)

print(diag(cor_var1_sup_cca_mcca) ^ 2)
print(diag(cor_var2_sup_cca_mcca) ^ 2)