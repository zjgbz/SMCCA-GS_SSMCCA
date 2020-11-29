rm(list=ls())
library(PMA)

args = commandArgs(trailingOnly=TRUE)

seed_list = c(1, 2, 3, 4, 5)
# seed_list = c(1)
vec_repeat_list = c(10)
matrix_num = 2
matrix_list = list()
noise_coef = as.numeric(args[1])
z_coef = as.numeric(args[2])
beta_coef = as.numeric(args[3])

df = data.frame(matrix(ncol = 0, nrow = 10))

nrow = 40
ncol_rel = 5
rep = 10
noise_sd = 10
ncol = ncol_rel
seed = 1
set.seed(seed)
z_matrix = matrix(rnorm(nrow * ncol, mean = 10, sd = 5), nrow, ncol_rel)
x1 = matrix(rnorm(nrow * ncol, mean = 3, sd = 2), nrow, ncol)
x2 = matrix(rnorm(nrow * ncol, mean = 3, sd = 2), nrow, ncol)

for (i in 1:ncol_rel) {
  z_i = z_matrix[, i]
  print(length(z_i))
  for (j in (i - 1) * rep + 1:i * rep) {
    a = rnorm(nrow * 1, mean = 0, sd = noise_sd)
    print(length(a))
    x1[, j] = z_i + rnorm(nrow * 1, mean = 0, sd = noise_sd)
    x2[, j] = z_i + rnorm(nrow * 1, mean = 0, sd = noise_sd)
  }
}

z_coef4y = rnorm(ncol_rel * 1, mean = 0, sd = 1)
z_combin = sweep(z_matrix, 2, z_coef4y, "*")
Y = rowSums(z_combin)

# moData = list(X1, X2)
outcome = "quantitative"
ncomp = 15
seed_list = c(1:15)

X1 = scale(X1_raw, T, T)
X2 = scale(X2_raw, T, T)

outcome = "quantitative"
num_components = 10

# built-in sup CCA
perm.out = CCA.permute(x=X1, z=X2, typex="standard",typez="standard",
                       outcome=outcome, y=y)
out <- CCA(x=X1,z=X2,typex="standard",typez="standard",K=num_components,
           penaltyx=perm.out$bestpenaltyx,penaltyz=perm.out$bestpenaltyz,
           v=perm.out$v.init, outcome=outcome, y=y)

w1_sup_cca = out$u
w2_sup_cca = out$v

var1_sup_cca = as.data.frame(as.matrix(X1) %*% w1_sup_cca)
var2_sup_cca = as.data.frame(as.matrix(X2) %*% w2_sup_cca)

# customized sup multi-CCA
source("sup_MultiCCA.R")

moData = list(X1, X2)

fullMCCA_1 = sup_MultiCCA.permute(moData, y, outcome)
fullMCCA2_1 = sup_MultiCCA(moData, y, outcome, penalty = fullMCCA_1$bestpenalties, ncomponents = num_components)
X1_sel = fullMCCA2_1$xlist[[1]]
X2_sel = fullMCCA2_1$xlist[[2]]
X1_sel_norm = scale(X1_sel, T, T)
X2_sel_norm = scale(X2_sel, T, T)

w1_sup_mcca = fullMCCA2_1$ws[[1]]
w2_sup_mcca = fullMCCA2_1$ws[[2]]

var1_sup_mcca = as.data.frame(as.matrix(X1_sel_norm) %*% w1_sup_mcca)
var2_sup_mcca = as.data.frame(as.matrix(X2_sel_norm) %*% w2_sup_mcca)

# compare the results
cor_var1_sup_cca_mcca = cor(var1_sup_cca, var1_sup_mcca)
cor_var2_sup_cca_mcca = cor(var2_sup_cca, var2_sup_mcca)

print("\n")
print(diag(cor_var1_sup_cca_mcca) ^ 2)
print(diag(cor_var2_sup_cca_mcca) ^ 2)
cor_var_df = cbind(as.data.frame(diag(cor_var1_sup_cca_mcca) ^ 2), diag(cor_var2_sup_cca_mcca) ^ 2)
colnames(cor_var_df) = c(sprintf("var1_%d", seed), sprintf("var2_%d", seed))
df = cbind(df, cor_var_df)
# built-in unsup multi-CCA

cor_var_filename = sprintf("noise_coef-%d_z_coef-%d_beta_coef-%d.tsv", noise_coef, z_coef, beta_coef)
cor_var_dir_filename = file.path(".", cor_var_filename)
write.table(df, file = cor_var_filename, row.names = FALSE, col.names = TRUE, sep = '\t')
