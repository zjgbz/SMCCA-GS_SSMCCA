# Functions below are developed by Min-Zhi Jiang

# Unsupervised Sparse Multiple Canonical Correlation Analysis
# with Gram–Schmidt Process

MultiCCA_GS <- function(moData, update_type="nores", opt_num=4, ncomponents=1, nperms=10, niter=25, cca_seed=42) {

	if (opt_num == "full" || opt_num > ncomponents) {
		opt_num = ncomponents
	}

	mat_num = length(moData)
	moData_scale = list()
	weight_list = list()
	cv_list = list()
	set.seed(cca_seed)
	for (mat_idx in 1:mat_num) {
		mat_i = moData[[mat_idx]]
		mat_i_scale = scale(mat_i, T, T)
		moData_scale[[mat_idx]] = mat_i_scale
		weight_list[[mat_idx]] = matrix(0, dim(mat_i)[2], ncomponents)
		cv_list[[mat_idx]] = matrix(0, dim(mat_i)[1], ncomponents)
	}

	penalty_list = list()
	niter_perm = 3
	
	for (comp_i in 1:ncomponents) {
		print(comp_i)
		if (comp_i <= opt_num) {
			fullMCCA = MultiCCA.permute(moData_scale, update_type=update_type, nperms=nperms, niter=niter_perm)
		}
		fullMCCA2 = MultiCCA(moData_scale, update_type=update_type, penalty=fullMCCA$bestpenalties, niter=niter, ncomponents=1)

		penalty_list[[comp_i]] = as.data.frame(fullMCCA$bestpenalties)
		# View(fullMCCA$bestpenalties)
		colnames(penalty_list[[comp_i]]) = comp_i

		moData_scale_old = moData_scale
		rm(moData_scale)
		moData_scale = list()
		for (mat_idx in 1:mat_num) {
			mat_i_scale_old = moData_scale_old[[mat_idx]]

			w_i = fullMCCA2$ws[[mat_idx]]
			var_i = mat_i_scale_old %*% w_i
			mat_i_scale = mat_i_scale_old - var_i %*% t(w_i)

			moData_scale[[mat_idx]] = mat_i_scale
			weight_list[[mat_idx]][, comp_i] = w_i
			cv_list[[mat_idx]][, comp_i] = var_i
		}
		print(sprintf("%d/%d canonical variables completed.", comp_i, ncomponents))
	}

	penalty_df = t(do.call(cbind, penalty_list))

	out = list(canon_var=cv_list, weight=weight_list, penalty=penalty_df)
	class(out) = "MultiCCA_GS"
	return (out)
}

# Supervised Sparse Multiple Canonical Correlation Analysis
# based on sup_MultiCCA_GS

sup_MultiCCA_GS <- function(moData, y, outcome, opt_num=4, update_type="nores", qt_list=NULL, ncomponents=1, nperms=10, niter=25, cca_seed=42) {
	# moData should be the list of dataframe rather than R matrix
	mat_num = length(moData)

	if (is.null(qt_list)) {
		qt_list = rep(0.8, mat_num)
	} else {
		if (length(qt_list) != mat_num) {
			stop("qt_list should include all threshold for each assay.")
		}
	}

	if (opt_num == "full" || opt_num > ncomponents) {
		opt_num = ncomponents
	}
	
	filter_out = MultiCCA.Phenotype.ZeroSome(moData, y, qt_list, cens=NULL, outcome=outcome)
	
	feature_dropped = filter_out$feature_dropped
	feature_kept = filter_out$feature_kept
	moData_scale_df = filter_out$xlist_sel
	moData_scale = df_list2matrix(moData_scale_df, mat_num)
	
	weight_sel_list = list()
	weight_list = list()
	cv_list = list()
	set.seed(cca_seed)
	for (mat_idx in 1:mat_num) {
		mat_i = moData_scale[[mat_idx]]
		mat_raw_i = moData[[mat_idx]]

		tmp_w_sel_mat = matrix(0, dim(mat_i)[2], ncomponents)
		tmp_w_sel = as.data.frame(tmp_w_sel_mat)
		rownames(tmp_w_sel) = colnames(mat_i)
		weight_sel_list[[mat_idx]] = tmp_w_sel

		tmp_w_mat = matrix(0, dim(mat_raw_i)[2], ncomponents)
		tmp_w = as.data.frame(tmp_w_mat)
		rownames(tmp_w) = colnames(mat_raw_i)
		weight_list[[mat_idx]] = tmp_w

		cv_list[[mat_idx]] = matrix(0, dim(mat_i)[1], ncomponents)
	}
	penalty_list = list()
	niter_perm=3
	
	for (comp_i in 1:ncomponents) {
		print(comp_i)
		if (comp_i <= opt_num) {
			fullMCCA = MultiCCA.permute(moData_scale, update_type=update_type, nperms=nperms, niter=niter_perm)
		}
		fullMCCA2 = MultiCCA(moData_scale, update_type=update_type, penalty=fullMCCA$bestpenalties, niter=niter, ncomponents=1)

		penalty_list[[comp_i]] = data.frame(fullMCCA$bestpenalties)
		colnames(penalty_list[[comp_i]]) = comp_i

		moData_scale_old = moData_scale
		rm(moData_scale)
		moData_scale = list()
		for (mat_idx in 1:mat_num) {
			mat_i_scale_old = moData_scale_old[[mat_idx]]

			w_i = fullMCCA2$ws[[mat_idx]]
			var_i = mat_i_scale_old %*% w_i
			mat_i_scale = mat_i_scale_old - var_i %*% t(w_i)

			moData_scale[[mat_idx]] = mat_i_scale
			rownames(w_i) = colnames(moData_scale[[mat_idx]])
			weight_sel_list[[mat_idx]][, comp_i] = w_i[rownames(weight_sel_list[[mat_idx]]), ]
			cv_list[[mat_idx]][, comp_i] = var_i

			weight_dropped = as.data.frame(matrix(0, nrow=length(feature_dropped[[mat_idx]]), ncol=1))
			rownames(weight_dropped) = feature_dropped[[mat_idx]]
			
			tmp_w_rbind = rbind(w_i, weight_dropped)
			weight_list[[mat_idx]][, comp_i] = tmp_w_rbind[rownames(weight_list[[mat_idx]]), ]
		}
		print(sprintf("%d/%d canonical variables completed.", comp_i, ncomponents))
	}

	penalty_df = t(do.call(cbind, penalty_list))

	weight_list_df = weight_list
	weight_sel_list_df = weight_sel_list
	rm(weight_list)
	rm(weight_sel_list)
	weight_list = list()
	weight_sel_list = list()
	for (mat_idx in 1:mat_num) {
		weight_list[[mat_idx]] = as.matrix(weight_list_df[[mat_idx]])
		weight_sel_list[[mat_idx]] = as.matrix(weight_sel_list_df[[mat_idx]])
	}

	# print(pvals)
	out = list(canon_var=cv_list, weight=weight_list, weight_select=weight_sel_list, feature_dropped=feature_dropped, penalty=penalty_df)
	class(out) = "sup_MultiCCA_GS"
	return (out)
}

# Function df_list2matrix is used for convert the
# input list of objects to dataframes

df_list2matrix <- function(xlist_input, K) {
	xlist = list()
	for (k in 1:K) {
		if (is.data.frame(xlist_input[[k]])) {
			xlist[[k]] = as.matrix(xlist_input[[k]])
		} else {
			xlist[[k]] = xlist_input[[k]]
		}
	}
	return (xlist)
}

# Functions below modifed from Witten and Tibshirani (2009) PMA R package
MultiCCA.Phenotype.ZeroSome <- function(xlist, y_raw, qt_list, cens=NULL, outcome=c("quantitative", "survival", "multiclass"), standarize=TRUE){
	outcome <- match.arg(outcome)
	K = length(xlist)
	assay_name_list=c(1:K)

	xlist_sel = list()
	score_list = list()
	feature_dropped = list()
	feature_kept = list()
	if (is.data.frame(y_raw)) {
		y_raw = y_raw[, 1]
	}

	if (is.list(y_raw) == FALSE) {
		ylist = list()
		for (k in 1:K) {
			ylist[[k]] = y_raw
		}
	} else {
		ylist = y_raw
	}
	for(k in 1:K) {
		tmp_x = xlist[[k]]
		qt = qt_list[k]
		y = ylist[[k]]
		if (standarize == TRUE) {
			tmp_x = scale(tmp_x, T, T)
		}
		if (outcome=="quantitative") {
			score.x <- quantitative.func(t(tmp_x)[,!is.na(y)],y[!is.na(y)])$tt
		} else if (outcome=="survival") {
			score.x <- cox.func(t(tmp_x)[,!is.na(y)],y[!is.na(y)],cens[!is.na(y)])$tt
		} else if (outcome=="multiclass") {
			score.x <- multiclass.func(t(tmp_x)[,!is.na(y)],y[!is.na(y)])$tt
		}

		keep.x <- abs(score.x) >= quantile(abs(score.x), qt)

		xnew <- tmp_x
		xlist_sel[[k]] = xnew[, keep.x]
		xlist_drop = xnew[,!keep.x]
		feature_dropped[[k]] = colnames(xlist_drop)
		feature_kept[[assay_name_list[k]]] = keep.x
		score_list[[k]] = score.x
	}
	
	return(list(xlist_sel = xlist_sel, feature_dropped = feature_dropped, feature_kept = feature_kept))
}

UpdateW <- function(xlist, i, K, sumabsthis, ws, ws.final, update_type="nores"){
	tots <- 0
	for(j in (1:K)[-i]){
		diagmat <- (t(ws.final[[i]])%*%t(xlist[[i]]))%*%(xlist[[j]]%*%ws.final[[j]])
		diagmat[row(diagmat)!=col(diagmat)] <- 0
		if (update_type == "res") {
			tots <- tots + t(xlist[[i]])%*%(xlist[[j]]%*%ws[[j]]) - ws.final[[i]]%*%(diagmat%*%(t(ws.final[[j]])%*%ws[[j]]))
		} else if (update_type == "nores") {
			tots <- tots + t(xlist[[i]])%*%(xlist[[j]]%*%ws[[j]])
		}
	}
	sumabsthis <- BinarySearch(tots, sumabsthis)
	w <- soft(tots, sumabsthis)/l2n(soft(tots, sumabsthis))
	return(w)
}

MultiCCA.permute <- function(xlist, update_type="nores", penalties=NULL, ws=NULL, nperms=10, niter=3, trace=TRUE, standardize=TRUE){
	call <- match.call()
	K <- length(xlist)
	for(k in 1:K){
		if(ncol(xlist[[k]])<2) stop("Need at least 2 features in each data set!")
		if(standardize) xlist[[k]] <- scale(xlist[[k]], T, T)
	}
	
	if(is.null(penalties)){
		penalties <- matrix(NA, nrow=K, ncol=10)
		for(k in 1:K){
			penalties[k,] <- pmax(seq(.1, .8, len=10)*sqrt(ncol(xlist[[k]])),1.1)
		}
	}
	numnonzeros <- NULL
	if(!is.matrix(penalties)) penalties <- matrix(1,nrow=K,ncol=1)%*%matrix(penalties,nrow=1)
	permcors <- matrix(NA, nrow=nperms, ncol=ncol(penalties))
	cors <- numeric(ncol(penalties))
	for(i in 1:ncol(penalties)){
		out <- MultiCCA(xlist, update_type=update_type, penalty=penalties[,i], niter=niter, ws=ws, trace=trace)
		cors[i] <- GetCors(xlist, out$ws, K)
		numnonzeros <- c(numnonzeros, sum(out$numnonzeros))
		ws.init  <- out$ws.init
	}
	cat(fill=TRUE)
	for(j in 1:nperms){
		if(trace) cat("Permutation ", j, "of " , nperms ,fill=TRUE)
		xlistperm <- xlist
		for(k in 1:K){
			xlistperm[[k]] <- xlistperm[[k]][sample(1:nrow(xlistperm[[k]])),]
		}
		for(i in 1:ncol(penalties)){
			out <- MultiCCA(xlistperm, update_type=update_type, penalty=penalties[,i], niter=niter, ws=ws, trace=FALSE)
			permcors[j,i] <- GetCors(xlistperm, out$ws, K)
		}
	}
	pvals =zs =  NULL
	for(i in 1:ncol(penalties)){
		pvals <- c(pvals, mean(permcors[,i]>=cors[i]))
		zs <- c(zs, (cors[i]-mean(permcors[,i]))/(sd(permcors[,i])+.05))
	}
	if(trace) cat(fill=TRUE)
	out <- list(pvals=pvals, zstat=zs, bestpenalties=penalties[,which.max(zs)], cors=cors, corperms=permcors, numnonzeros=numnonzeros, ws.init=ws.init, call=call, penalties=penalties, nperms=nperms)
	class(out) <- "MultiCCA.permute"
	return(out)
}

MultiCCA <- function(xlist, update_type, penalty=NULL, ws=NULL, niter=25, ncomponents=1, trace=TRUE, standardize=TRUE){
	for(i in 1:length(xlist)){
		if(ncol(xlist[[i]])<2) stop("Need at least 2 features in each data set.")
	}
	call <- match.call()
	K <- length(xlist)

	for(k in 1:K){
		if(standardize) xlist[[k]] <- scale(xlist[[k]], T, T)
	}
	if(!is.null(ws)){
		makenull <- FALSE
		for(i in 1:K){
			if(ncol(ws[[i]])<ncomponents) makenull <- TRUE
		}
		if(makenull) ws <- NULL
	}
	if(is.null(ws)){
		ws <- list()
		for(i in 1:K) ws[[i]] <- matrix(svd(xlist[[i]])$v[,1:ncomponents], ncol=ncomponents)
	}
	if(is.null(penalty)){
		penalty <- rep(4, K) # this is the default value of sumabs
	}
	ws.init <- ws
	if(length(penalty)==1) penalty <- rep(penalty, K)
	
	for(i in 1:length(xlist)){
		if(penalty[i]>sqrt(ncol(xlist[[i]]))) stop("L1 bound of weights should be no more than sqrt of the number of columns of the corresponding data set.", fill=TRUE)
	}
	ws.final <- list()
	for(i in 1:length(ws)) ws.final[[i]] <- matrix(0, nrow=ncol(xlist[[i]]), ncol=ncomponents)
	cors <- NULL
	for(comp in 1:ncomponents){
		ws <- list()
		for(i in 1:length(ws.init)) ws[[i]] <- ws.init[[i]][,comp]
		curiter <- 1
		crit.old <- -10
		crit <- -20
		storecrits <- NULL
		while(curiter<=niter && abs(crit.old-crit)/abs(crit.old)>.001 && crit.old!=0){
			crit.old <- crit
			crit <- GetCrit(xlist, ws, K)
			storecrits <- c(storecrits,crit)
			if(trace) cat(curiter, fill=FALSE)
			curiter <- curiter+1
			for(i in 1:K){
				ws[[i]] <- UpdateW(xlist, i, K, penalty[i], ws, ws.final, update_type)
			}
		}
		for(i in 1:length(ws)) ws.final[[i]][,comp] <- ws[[i]]
		cors <- c(cors, GetCors(xlist, ws,K))
	}
	out <- list(ws=ws.final, ws.init=ws.init, K=K, call=call, penalty=penalty, cors=cors)
	class(out) <- "MultiCCA"
	return(out)
}

# Functions below inherits from Witten and Tibshirani (2009) PMA R package

soft <- function(x,d){
	return(sign(x)*pmax(0, abs(x)-d))
}

mean.na <- function(vec){
	return(mean(vec[!is.na(vec)]))
}

l2n <- function(vec){
	a <- sqrt(sum(vec^2))
	if(a==0) a <- .05
	return(a)
}

GetCrit <- function(xlist, ws, K){
	crit <- 0
	for(i in 2:K){
		for(j in 1:(i-1)){
				crit <- crit + t(ws[[i]])%*%t(xlist[[i]])%*%xlist[[j]]%*%ws[[j]]
		}
	}
	return(crit)
}

GetCors <- function(xlist, ws, K){
	cors <- 0
	for(i in 2:K){
		for(j in 1:(i-1)){
			thiscor  <-  cor(xlist[[i]]%*%ws[[i]], xlist[[j]]%*%ws[[j]])
			if(is.na(thiscor)) thiscor <- 0
			cors <- cors + thiscor
		}
	}
	return(cors)
}

ftrans <- function(x){ return(.5*log((1+x)/(1-x))) }

BinarySearch <- function(argu,sumabs){
	if(l2n(argu)==0 || sum(abs(argu/l2n(argu)))<=sumabs) return(0)
	lam1 <- 0
	lam2 <- max(abs(argu))-1e-5
	iter <- 1
	while(iter < 150){
		su <- soft(argu,(lam1+lam2)/2)
		if(sum(abs(su/l2n(su)))<sumabs){
			lam2 <- (lam1+lam2)/2
		} else {
			lam1 <- (lam1+lam2)/2
		}
		if((lam2-lam1)<1e-6) return((lam1+lam2)/2)
		iter <- iter+1
	}
	warning("Didn't quite converge")
	return((lam1+lam2)/2)
}

varr <- function(x, meanx=NULL){
	n <- ncol(x)
	p <- nrow(x)
	Y <-matrix(1,nrow=n,ncol=1)
	if(is.null(meanx)){   meanx <- rowMeans(x)}
	ans<- rep(1, p)
	xdif <- x - meanx %*% t(Y)
	ans <- (xdif^2) %*% rep(1/(n - 1), n)
	ans <- drop(ans)
	return(ans)
}

quantitative.func  <- function(x,y,s0=0){
	
	# regression of x on y

	my=mean(y)
	yy <- y-my
	temp <- x%*%yy
	mx=rowMeans(x)
	syy= sum(yy^2)

	scor <- temp/syy
	b0hat <- mx-scor*my
	ym=matrix(y,nrow=nrow(x),ncol=ncol(x),byrow=T)
	xhat <- matrix(b0hat,nrow=nrow(x),ncol=ncol(x))+ym*matrix(scor,nrow=nrow(x),ncol=ncol(x))
	sigma <- sqrt(rowSums((x-xhat)^2)/(ncol(xhat)-2))
	sd <- sigma/sqrt(syy)
	tt <- scor/(sd+s0)

	return(list(tt=tt, numer=scor, sd=sd))
}

multiclass.func <- function(x,y,s0=0){

	##assumes y is coded 1,2...

	nn <- table(y)
	m <- matrix(0,nrow=nrow(x),ncol=length(nn))
	v <- m
	for(j in 1:length(nn)){
		m[,j] <- rowMeans(x[,y==j])
		v[,j] <- (nn[j]-1)*varr(x[,y==j], meanx=m[,j])
	}
	mbar <- rowMeans(x)
	mm <- m-matrix(mbar,nrow=length(mbar),ncol=length(nn))
	fac <- (sum(nn)/prod(nn))
	scor <- sqrt(fac*(apply(matrix(nn,nrow=nrow(m),ncol=ncol(m),byrow=TRUE)*mm*mm,1,sum)))

	sd <- sqrt(rowSums(v)*(1/sum(nn-1))*sum(1/nn))
	tt <- scor/(sd+s0)
	mm.stand=t(scale(t(mm),center=FALSE,scale=sd))
	return(list(tt=tt, numer=scor, sd=sd,stand.contrasts=mm.stand))
}
