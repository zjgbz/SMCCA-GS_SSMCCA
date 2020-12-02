rm(list = ls())
start_time <- Sys.time()
# library(PMA)
library(ggplot2)
library(tidyr)
library(reshape2)
library(corrplot)
source("sup_MultiCCA.R")

get_cor <- function(xlist_input, ws, K){
	xlist = df_list2matrix(xlist_input, K)
	cors <- 0
	for(i in 2:K){
		for(j in 1:(i-1)){
			thiscor <- cor(xlist[[i]]%*%ws[[i]], xlist[[j]]%*%ws[[j]])
			if (dim(ws[[i]])[2] > 1) {
				thiscor_array = diag(thiscor)
			} else {
				thiscor_array = thiscor
			}
			# if(is.na(thiscor_array)) thiscor_array <- 0
			cors <- cors + thiscor_array
		}
	}
	return(cors)
}

count_contrib <- function(wi_i, rep, ncol_rel) {
	wi_i = as.matrix(wi_i)
	contrib_idx = row(wi_i)[which(!wi_i == 0)]
	contrib_cum = data.frame(matrix(0, nrow = 1, ncol = ncol_rel + 1))
	for (rel_i in 1:ncol_rel) {
		contrib_cum[, rel_i] = sum(contrib_idx <= rel_i * rep & contrib_idx > (rel_i - 1) * rep)
		colnames(contrib_cum)[rel_i] = sprintf("group_%d", rel_i)
	}
	contrib_cum[, ncol_rel + 1] = sum(contrib_idx > rel_i * rep)
	colnames(contrib_cum)[ncol_rel + 1] = "noise"
	return (contrib_cum)
}

organize_contrib <- function(rep, ncol_rel, w) {
	contrib_w1 = apply(w1, 2, count_contrib, rep = rep, ncol_rel = ncol_rel)
	contrib_w1_df = do.call("rbind", contrib_w1)
	contrib_w1_forbar_df = t(contrib_w1_df)
	ncomp = dim(contrib_w1_forbar_df)[2]
	cv_list = c()
	for (cv_idx in 1:ncomp) {
		cv_list[cv_idx] = cv_idx
	}
	colnames(contrib_w1_forbar_df) = cv_list
	contrib_w1_forbar = as.matrix(contrib_w1_forbar_df)
}

save_cor_plot <- function(matrix4cor, tl_cex, cl_cex, output_dir, filename_suffix, filename_prefix = NULL) {
	if (is.null(filename_prefix)) {
		filename = filename_suffix
	} else {
		filename = sprintf("%s_%s.png", filename_prefix, filename_suffix)
	}
	cor_mat = cor(matrix4cor, matrix4cor)
	cor_mat_dir_filename = file.path(output_dir, filename)
	cor_mat_abs = abs(cor_mat)
	order.hc2 <- corrMatOrder(cor_mat_abs, order = "hclust", hclust.method = "ward.D")
	M.hc2 <- cor_mat[order.hc2,order.hc2]
	melted_cormat <- melt(M.hc2)
	ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + theme(
		axis.text.x = element_text(angle = 90, size = 12),
		axis.text.y = element_text(angle = 0, size = 12),
		axis.title.x = element_blank(),
		axis.title.y = element_blank()) +
	scale_fill_gradient2(low="blue", mid="lightblue", high="red") +
	geom_tile()
	ggsave(file=cor_mat_dir_filename)
}

assay_cor_plot <- function(matrix4cor, rep, output_dir, filename_suffix, filename_prefix = NULL) {
	if (is.null(filename_prefix)) {
		filename = filename_suffix
	} else {
		filename = sprintf("%s_%s.png", filename_prefix, filename_suffix)
	}
	cor_mat = cor(matrix4cor, matrix4cor)
	melted_cormat <- melt(cor_mat)
	cor_mat_dir_filename = file.path(output_dir, filename)
	# png(file = cor_mat_dir_filename)
	# col_mat <- colorRampPalette(c("blue", "white", "red"))(20)
	# heatmap(x = cor_mat, col = col_mat, symm = TRUE)
	# dev.off()

	ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + theme(
		axis.text.x = element_text(angle = 90, size = 12),
		axis.text.y = element_text(angle = 0, size = 12),
		axis.title.x = element_blank(),
		axis.title.y = element_blank()) +
	scale_x_continuous(breaks=seq(1, dim(cor_mat)[1], rep)) + scale_y_continuous(breaks=seq(1, dim(cor_mat)[1], rep)) +
	scale_fill_gradient2(low="blue", mid="lightblue", high="red") +
	geom_tile()
	ggsave(file=cor_mat_dir_filename)
}

cv_assay_cor_plot <- function(assay, cv, rep, output_dir, filename_suffix, filename_prefix = NULL) {
	if (is.null(filename_prefix)) {
		filename = filename_suffix
	} else {
		filename = sprintf("%s_%s.png", filename_prefix, filename_suffix)
	}
	cor_mat = cor(assay, cv)
	melted_cormat <- melt(cor_mat)
	cor_mat_dir_filename = file.path(output_dir, filename)
	ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + theme(
		axis.text.x = element_text(angle = 90, size = 12),
		axis.text.y = element_text(angle = 0, size = 12),
		axis.title.x = element_blank(),
		axis.title.y = element_blank()) +
	scale_x_continuous(breaks=seq(1, dim(x1)[2], rep)) +
	scale_fill_gradient2(low="blue", mid="lightblue", high="red") +
	geom_tile()
	ggsave(file=cor_mat_dir_filename, width=20, height=10)
}

save_weight <- function(w_mat, output_dir, filename_suffix, filename_prefix = NULL) {
	if (is.null(filename_prefix)) {
		filename = filename_suffix
	} else {
		filename = sprintf("%s_%s.tsv", filename_prefix, filename_suffix)
	}
	w_dir_filename = file.path(output_dir, filename)
	write.table(w_mat, file = w_dir_filename, row.names=FALSE, col.names=FALSE, sep='\t')
}

rehead_mat <- function(mat, head_prefix) {
	head_list = c()
	col_num = dim(mat)[2]
	for (col_i in 1:col_num) {
		head_list[col_i] = sprintf("%s_%d", head_prefix, col_i)
	}
	colnames(mat) = head_list
	return (mat)
}

# args = commandArgs(trailingOnly=TRUE)
# cca_seed = as.numeric(args[1])
# feature_factor = as.numeric(args[2])

# cca_seed = 1
feature_factor = 1

size_factor = 100
nrow = 4 * size_factor
ncol_rel = 5
rep = 1 * size_factor * feature_factor
noise_sd = 5
ncol = 8 * size_factor * feature_factor
# matrix_seed_vec = c(1)
cca_seed_vec = c(1:20)
cca_seed_vec = c(1)
noise_mode = "noise-reinf"
# noise_mode = "noise-keep"
cca_fashion = "sup-mcca"
ncomp = 15

# x1_feature_group_list = c("x1_g1", "x1_g2", "x1_g3", "x1_g4", "x1_g5", "x1_g6", "x1_g7", "x1_g8")
# x2_feature_group_list = c("x2_g1", "x2_g2", "x2_g3", "x2_g4", "x2_g5", "x2_g6", "x2_g7", "x2_g8")
# x1_x2_feature_group_list = c(x1_feature_group_list, x2_feature_group_list)

for (matrix_seed in matrix_seed_vec) {
	for (cca_seed in cca_seed_vec) {
		output_dir = file.path(".", sprintf("matrix-seed%d_cca-seed%d_size%d_dim%d_%s", matrix_seed, cca_seed, size_factor, feature_factor, noise_mode))
		dir.create(output_dir, showWarnings = TRUE, recursive = TRUE, mode = "0777")
		existed_image_prefix = sprintf("%s_cca-seed%d_matrix-seed%d_nrow%d_ncol-rel%d_rep%d_noise-sd%d_%s_ncol%d", cca_fashion, cca_seed, matrix_seed, nrow, ncol_rel, rep, noise_sd, noise_mode, ncol)
		existed_image_filename = sprintf("%s.RData", existed_image_prefix)
		existed_image_dir_filename = file.path(output_dir, existed_image_filename)

		if (!file.exists(existed_image_dir_filename)) {
			set.seed(matrix_seed)
			z_matrix = matrix(rnorm(nrow * ncol, mean = 10, sd = 5), nrow, ncol_rel)
			x1 = matrix(rnorm(nrow * ncol, mean = 3, sd = 2), nrow, ncol)
			x2 = matrix(rnorm(nrow * ncol, mean = 3, sd = 2), nrow, ncol)
			y_coef = matrix(runif(ncol_rel, min = 0, max = 1), ncol_rel, 1)
			y = z_matrix %*% y_coef

			for (i in 1:ncol_rel) {
				z_i = z_matrix[, i]
				start = (i - 1) * rep + 1
				end = i * rep
				if (noise_mode == "noise-reinf") {
					noise_assign = i * noise_sd / 2
				} else {
					noise_assign = noise_sd
				}
				for (j in start:end) {
					x1[, j] = z_i + rnorm(nrow * 1, mean = 0, sd = noise_assign)
					x2[, j] = z_i + rnorm(nrow * 1, mean = 0, sd = noise_assign)
				}
			}
			x1_df = as.data.frame(x1)
			x2_df = as.data.frame(x2)

			set.seed(cca_seed)

			input_list = list(x1_df, x2_df)
			# x1_x2 = cbind(x1, x2)
			# assay_cor_plot(x1, rep, output_dir, "cor-x1", existed_image_prefix)
			# assay_cor_plot(x2, rep, output_dir, "cor-x2", existed_image_prefix)
			# assay_cor_plot(x1_x2, rep, output_dir, "cor-x1-x2", existed_image_prefix)

			fullMCCA = sup_MultiCCA.permute(xlist_raw = input_list, y = y, outcome = "quantitative")
			input_list_select = fullMCCA$xlist
			save.image(file=existed_image_dir_filename)
		} else if (file.exists(existed_image_dir_filename)) {
			load(existed_image_dir_filename)
			print("image loaded")
		}

		K = length(input_list)

		nperm_test = 10
		permcors <- matrix(NA, nrow = nperm_test, ncol = ncomp)
		perm_seed = 1
		nrep = 1
		nperm = 10

		existed_perm_filename = sprintf("%s_ncomp%d_nperm%d_perm-seed%d.RData", existed_image_prefix, ncomp, nperm_test, perm_seed)
		existed_perm_dir_filename = file.path(output_dir, existed_perm_filename)
		if (!file.exists(existed_perm_dir_filename)) {
			cors = matrix(NA, nrow = 1, ncol = ncomp)
			out = sup_MultiCCA(input_list, input_list_select, fullMCCA$feature_dropped, penalty = fullMCCA$bestpenalties, ncomponents = ncomp)
			fullcca_time <- Sys.time()

			w1 = out$ws[[1]]
			w2 = out$ws[[2]]

			save_weight(w1, output_dir, "w1", filename_prefix = existed_image_prefix)
			save_weight(w2, output_dir, "w2", filename_prefix = existed_image_prefix)
			
			mycols = c("white","yellow", "orange", "red", "purple", "blue")
			col = rep(mycols, ncomp)

			contrib_w1_forbar_filename = sprintf("contrib_w1_%s_n%d_var-%s_noise%d_%s_ncomp%d_matrix-seed%d_cca-seed%d.png", cca_fashion, ncomp, length(input_list), noise_sd, noise_mode, ncomp, matrix_seed, cca_seed)
			contrib_w1_forbar_dir_filename = file.path(output_dir, contrib_w1_forbar_filename)
			png(contrib_w1_forbar_dir_filename)
			contrib_w1_forbar = organize_contrib(rep, ncol_rel, w1)
			barplot(contrib_w1_forbar, col = col)
			legend(x = "right", legend = rownames(contrib_w1_forbar), fill = mycols)
			dev.off()
			
			contrib_w2_forbar_filename = sprintf("contrib_w2_%s_n%d_var-%s_noise%d_%s_ncomp%d_matrix-seed%d_cca-seed%d.png", cca_fashion, ncomp, length(input_list), noise_sd, noise_mode, ncomp, matrix_seed, cca_seed)
			contrib_w2_forbar_dir_filename = file.path(output_dir, contrib_w2_forbar_filename)		
			png(contrib_w2_forbar_dir_filename)
			contrib_w2_forbar = organize_contrib(rep, ncol_rel, w2)
			barplot(contrib_w2_forbar, col = col)
			legend(x = "right", legend = rownames(contrib_w2_forbar), fill = mycols)
			dev.off()
	
			cv1_raw = scale(x1, T, T) %*% w1
			cv2_raw = scale(x2, T, T) %*% w2
			cv1 = rehead_mat(cv1_raw, "cv1")
			cv2 = rehead_mat(cv2_raw, "cv2")
			cv1_cv2 = cbind(cv1, cv2)
			save_cor_plot(cv1, tl_cex = 3, cl_cex = 1.5, output_dir, "cor-cv1", existed_image_prefix)
			save_cor_plot(cv2, tl_cex = 3, cl_cex = 1.5, output_dir, "cor-cv2", existed_image_prefix)
			save_cor_plot(cv1_cv2, tl_cex = 0.9, cl_cex = 1.5, output_dir, "cor-cv1-cv2", existed_image_prefix)

			cv_assay_cor_plot(x1, cv1, rep, output_dir, "cor-cv1-x1", existed_image_prefix)
			cv_assay_cor_plot(x2, cv2, rep, output_dir, "cor-cv2-x2", existed_image_prefix)

			cors[1, ] <- get_cor(input_list, out$ws, K)

			set.seed(perm_seed)
			for (perm_i in 1:nperm_test) {
			    xlistperm <- input_list
			    for(k in 1:K){
			      xlistperm[[k]] <- xlistperm[[k]][sample(1:nrow(xlistperm[[k]])),]
			    }
				out_perm = sup_MultiCCA(input_list, xlistperm, fullMCCA$feature_dropped, penalty = fullMCCA$bestpenalties, ncomponents = ncomp)
				permcors[perm_i, ] <- get_cor(xlistperm, out_perm$ws, K)
			}

			pvals = NULL
			for (i in 1:ncomp) {
				pvals <- c(pvals, mean(permcors[,i] >= cors[i]))
			}
			save.image(file=existed_perm_dir_filename)
		} else if (file.exists(existed_perm_dir_filename)) {
			load(existed_perm_dir_filename)
			print("pvalue image loaded")
		}
		print(pvals)

		pvals.m <- melt(pvals)
		ggplot(pvals.m, aes(x=seq_along(value), y = value)) + coord_cartesian(ylim = c(-0.05, 1.05)) + scale_x_discrete(limits=factor(c(1:ncomp))) +
			geom_line(aes()) + labs(y = "p-value", x = "#CV",title = "p-value vs #CV") + theme(
			title =element_text(size=30, face='bold'),
			axis.text.x = element_text(face="bold", size=20),
			axis.text.y = element_text(face="bold", size=20),
			axis.title.x = element_text(face="bold", size=20),
			axis.title.y = element_text(face="bold", size=20),
			legend.title=element_text(size=30), 
			legend.text=element_text(size=20)
			)

		fig1_filename = sprintf("p-value_%s_n%d_var-%s_noise%d_%s_rep-%d_nperm%d_nperm-test%d_ncomp%d_matrix-seed%d_cca-seed%d_perm-seed%d.png", cca_fashion, ncomp, length(input_list), noise_sd, noise_mode, nrep, nperm, nperm_test, ncomp, matrix_seed, cca_seed, perm_seed)
		fig1_dir_filename = file.path(output_dir, fig1_filename)
		ggsave(file=fig1_dir_filename, width=35, height=10)

		cors.m <- melt(cors)
		ggplot(cors.m, aes(x=seq_along(value), y = value)) + coord_cartesian(ylim = c(-1.05, 1.05)) + scale_x_discrete(limits=factor(c(1:ncomp))) +
			geom_line(aes()) + labs(y = "sum cor", x = "#CV",title = "sum correlation vs #CV") + theme(
			title =element_text(size=30, face='bold'),
			axis.text.x = element_text(face="bold", size=20),
			axis.text.y = element_text(face="bold", size=20),
			axis.title.x = element_text(face="bold", size=20),
			axis.title.y = element_text(face="bold", size=20),
			legend.title=element_text(size=30), 
			legend.text=element_text(size=20)
			)
		fig2_filename = sprintf("sum_cor_%s_n%d_var-%s_noise%d_%s_rep-%d_nperm%d_nperm-test%d_ncomp%d_matrix-seed%d_cca-seed%d_perm-seed%d.png", cca_fashion, ncomp, length(input_list), noise_sd, noise_mode, nrep, nperm, nperm_test, ncomp, matrix_seed, cca_seed, perm_seed)
		fig2_dir_filename = file.path(output_dir, fig2_filename)
		ggsave(file=fig2_dir_filename, width=35, height=10)		
	}
}
end_time <- Sys.time()
fullcca_running_time = fullcca_time - start_time
print(fullcca_running_time)
full_running_time = end_time - start_time
print(full_running_time)
