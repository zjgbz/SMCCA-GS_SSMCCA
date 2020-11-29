rm(list = ls())

library(PMA)
library(ggplot2)
library(tidyr)
library(reshape2)
library(dplyr)
library(plotly)
library(gridExtra)
library(nnet)
source("sup_MultiCCA.R")

store_result <- function(input_list, cca_fashion, ncomp, seed_list, nperm, y = NULL, outcome = NULL) {
  result1_list = list()
  result2_list = list()
  seed_num = length(seed_list)
  for (idx in 1:seed_num) {
    seed = seed_list[idx]
    set.seed(seed)
    if (cca_fashion == "scca") {
      X1 = input_list[[1]]
      X2 = input_list[[2]]
      perm.out = CCA.permute(x = X1, z = X2, typex = "standard", typez = "standard", nperms = nperm,
                             outcome = outcome, y = y)
      out <- CCA(x=X1, z=X2,typex="standard",typez="standard",K = ncomp,
                 penaltyx=perm.out$bestpenaltyx,penaltyz=perm.out$bestpenaltyz,
                 v=perm.out$v.init, outcome=outcome, y=y)
      w1 = out$u
      w2 = out$v
    } else if (cca_fashion == "unsup_mcca") {
      fullMCCA = MultiCCA.permute(input_list, nperms = nperm)
      fullMCCA2 = MultiCCA(input_list, penalty = fullMCCA$bestpenalties, ncomponents = ncomp)
      w1 = fullMCCA2$ws[[1]]
      w2 = fullMCCA2$ws[[2]]
    } else if (cca_fashion == "sup_mcca") {
      fullMCCA = sup_MultiCCA.permute(input_list, y, outcome, nperms = nperm)
      fullMCCA2 = sup_MultiCCA(input_list, y, outcome, penalty = fullMCCA$bestpenalties, ncomponents = ncomp)
      w1 = fullMCCA2$ws[[1]]
      w2 = fullMCCA2$ws[[2]]
    }
    norm1 = scale(input_list[[1]], T, T)
    norm2 = scale(input_list[[2]], T, T)
    var1 = as.data.frame(as.matrix(norm1) %*% w1)
    var2 = as.data.frame(as.matrix(norm2) %*% w2)
    result1_list[[idx]] = var1
    result2_list[[idx]] = var2
  }
  result_list = list(result1_list, result2_list)
  return (result_list)
}

generate_colname <- function(cor_df) {
  col_num = dim(cor_df)[2]
  for (col_i in 1:col_num) {
    colnames(cor_df)[col_i] = sprintf("CV_%02d", col_i)
  }
  return (cor_df)
}

iter_result <- function(seed_list, result_list, ncomp) {
  seed_num = length(seed_list)
  iter_idx_df = t(combn(seed_num, 2))
  # print(iter_idx_df)
  iter_num = dim(iter_idx_df)[1]
  result_num = length(result_list)
  cor_df_list = list()
  for (result_i in 1:result_num) {
    result = result_list[[result_i]]
    cor_df = data.frame(matrix(nrow = 0, ncol = ncomp))
    for (iter_i in 1:iter_num) {
      var_1_idx = iter_idx_df[iter_i, 1]
      var_2_idx = iter_idx_df[iter_i, 2]
      var_1 = result[[var_1_idx]]
      var_2 = result[[var_2_idx]]
      cor_var_1_2 = diag(cor(var_1, var_2)) ^ 2
      cor_df = rbind(cor_df, cor_var_1_2)
      cor_df_colname = generate_colname(cor_df)
    }
    cor_df_list[[result_i]] = cor_df_colname
  }
  return (cor_df_list)
}

violin_plot <- function(df, idx, seed_num, nperm, noise_sd) {
  df_violin_plot = df %>%
  gather(key="CV", value="rsquare") %>%
  ggplot(aes(x = CV, y = rsquare, fill = CV), show.legend = FALSE) + geom_violin() + labs(title = sprintf("var %d, %d rep, %d noise, %d nperm", idx, seed_num, noise_sd, nperm)) + theme(
    axis.text.x = element_text(face="bold", size=30),
    axis.text.y = element_text(face="bold", size=20),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    plot.title = element_text(size=22)
    )
  return (df_violin_plot)
}

violin_plot_list <- function(df_list, ncomp, df_num, seed_num, nperm, noise_sd) {
  df_num = length(df_list)
  df_violin_plot_list = list()
  for (idx in 1:df_num) {
    df = df_list[[idx]]
    df_violin_plot = violin_plot(df, idx, seed_num, nperm, noise_sd)
    df_violin_plot_list[[idx]] = df_violin_plot
  }
  fig_num = df_num
  fig_col = grid.arrange(grobs = df_violin_plot_list, layout_matrix = rbind(1,df_num), ncol=1, nrow=df_num)
  ggsave(file=sprintf("stability_%s_n%d_var-%s_noise%d_seed-%d_nperm%d.png", cca_fashion, ncomp, df_num, noise_sd, seed_num, nperm), plot=fig_col, width=35, height=10)
  dev.off()
  return (df_violin_plot_list)
}

sum_cor <- function(df_list) {
  ncomp = dim(df_list[[1]])[2]
  df_num = length(df_list)
  iter_idx_df = t(combn(df_num, 2))
  iter_num = dim(iter_idx_df)[1]
  sum_correlation_vec = rep(0, iter_num)
  for (iter_i in 1:iter_num){
    iter_idx_1 = iter_idx_df[iter_i, 1]
    iter_idx_2 = iter_idx_df[iter_i, 2]
    df1 = df_list[[iter_idx_1]]
    df2 = df_list[[iter_idx_2]]
    cor_df1_df2 = cor(df1, df2)
    tmp_cor = diag(cor_df1_df2)
    sum_correlation_vec =  sum_correlation_vec + tmp_cor
  }
  return (sum_correlation_vec)
}

sum_cor_list <- function(df_list_list) {
  rep_num = length(df_list_list)
  ncomp = dim(df_list_list[[1]][[1]])[2]
  sum_cor_df = matrix(NA, rep_num, ncomp)
  for (rep_i in 1:rep_num) {
    df_list = df_list_list[[rep_i]]
    sum_cor_df[rep_i, ] = sum_cor(df_list)
  }
  return (sum_cor_df)
}

reorg_result <- function(result_list) {
  var_num = length(result_list)
  rep_num = length(result_list[[1]])
  reorg_result_list = list()
  for (rep_i in 1:rep_num) {
    var_list = list()
    for (var_i in 1:var_num) {
      var_list[[var_i]] = result_list[[var_i]][[rep_i]]
    }
    reorg_result_list[[rep_i]] = var_list
  }
  return (reorg_result_list)
}

nrow = 40
ncol_rel = 5
rep = 10
noise_sd = 1
ncol = ncol_rel * rep + ncol_rel * rep
seed = 1
set.seed(seed)
z_matrix = matrix(rnorm(nrow * ncol, mean = 10, sd = 5), nrow, ncol_rel)
x1 = matrix(rnorm(nrow * ncol, mean = 3, sd = 2), nrow, ncol)
x2 = matrix(rnorm(nrow * ncol, mean = 3, sd = 2), nrow, ncol)

for (i in 1:ncol_rel) {
  z_i = z_matrix[, i]
  start = (i - 1) * rep + 1
  end = i * rep
  noise_assign = i * noise_sd
  for (j in start:end) {
    x1[, j] = z_i + rnorm(nrow * 1, mean = 0, sd = noise_assign)
    x2[, j] = z_i + rnorm(nrow * 1, mean = 0, sd = noise_assign)
  }
}

z_coef4y = rnorm(ncol_rel * 1, mean = 0, sd = 1)
z_combin = sweep(z_matrix, 2, z_coef4y, "*")
Y = rowSums(z_combin)

# moData = list(X1, X2)
outcome = "quantitative"
ncomp = 15
# seed_list = c(1:5)
seed_list = c(1:15)

# m1_col = 1000
# m2_col = 1000
# row_num = 400
# set.seed(100)
# m1 = matrix(rnorm(row_num*m1_col),row_num,m1_col) 
# m2 = matrix(rexp(row_num*m2_col, rate=.1), row_num,m2_col)
input_list = list(x1, x2)

cca_fashion = "scca"

# nperm_list = c(10, 50, 100)
nperm_list = c(10)
for (nperm in nperm_list) {
  result_list = store_result(input_list, cca_fashion, ncomp, seed_list, nperm, y = Y, outcome = outcome)
  reorg_result_list = reorg_result(result_list)
  sum_cor_df = sum_cor_list(reorg_result_list)
  
  df_melted_raw = melt(sum_cor_df)
  df_melted = data.frame(rep = df_melted_raw[["Var1"]], cv = df_melted_raw[["Var2"]], sum_corsq = df_melted_raw[["value"]])
  ggplot(df_melted, aes(x = cv, y = sum_corsq, group = rep, color = rep)) + geom_line(aes(color = 'cat')) + theme(
      axis.text.x = element_text(face="bold", size=30),
      axis.text.y = element_text(face="bold", size=20),
      axis.title.x = element_text(face="bold", size=30),
      axis.title.y = element_text(face="bold", size=20),
      legend.title=element_text(size=30), 
      legend.text=element_text(size=20)
      )
  ggsave(file=sprintf("sum_cor_%s_n%d_var-%s_noise%d_rep-%d_nperm%d.png", cca_fashion, ncomp, length(input_list), noise_sd, length(seed_list), nperm), width=35, height=10)
  # while (!is.null(dev.list())) dev.off()
  # dev.off()
  
  iter_cor = iter_result(seed_list, result_list, ncomp)
  violin_plot_list(iter_cor, ncomp, length(input_list), length(seed_list), nperm, noise_sd)
}
