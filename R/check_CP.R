#' Reconstruct the matching based on the estimated parameters
#' 
#' Calculate the joint density of the observed features of the pairs between agents
#' using the estimated parameters, given the feature matrices of the two sides. The 
#' observed joint density is also computed.  
#' 
#' @param ff formula; an \code{\link{formula}} object, of the form \code{
#' ~ <model terms>}. For the details on the possible \code{<model terms>}, see
#' \code{\link{rpm-terms}}.
#' @param theta The estimated parameters.
#' @param mu The observed matching matrix, where 1 represent a pairing, 0 otherwise. 
#' Each row is a woman, each column is a man. The order of the rows (columns)
#' needs to be the same as in \code{X} (\code{Z}). 
#' @param X Feature matrix for women. Each row is a woman, each column is a feature. 
#' The number of column is assumed to be the same as in \code{Zdata}. 
#' @param Z Feature matrix for men. Each row is a man, each column is a feature. 
#' The number of column is assumed to be the same as in \code{Xdata}.
#' @param symmetric If set to \code{TRUE}, the same utility coefficients are to be used 
#' for both the women's and men's sides; otherwise, separate estimates need to be provided for each side. 
#' @return This function returns a list consisting of the following elements: 
#' \item{pmfj_est}{The reconstructed joint density matrix based on the estimated parameters.}
#' \item{pmfj_obs}{The observed joint density matrix based on the data \code{X} and \code{Z}.}
#' @seealso fitrpm_R_CP
#' @references Menzel, Konrad (2015).
#' \emph{Large Matching Markets as Two-Sided Demand Systems}
#' Econometrica, Vol. 83, No. 3 (May, 2015), 897-941.
#' @keywords models
#' @examples
#' 
check_CP_latent = function(ff, theta, mu, X, Z, symmetric){
  
  n=nrow(X)+nrow(Z)
  
  Xdata <- cbind(1, X)
  Zdata <- cbind(1, Z)
  colnames(Xdata)[1] <- "Int"
  colnames(Zdata)[1] <- "Int"
  
  model.terms <- rownames(attr(terms.formula(ff), "factors"))
  temp <- strsplit(model.terms, "[(]")
  model.terms.names <- unlist(lapply(temp, `[[`, 1))
  temp2 <- gsub(")", "", gsub("\"", "", unlist(temp)))
  temp2 <- gsub("[ \t\n\r\f\v]", "", temp2)
  model.terms.coef.names <- temp2[-which(temp2 %in% model.terms.names)]
  model.terms.coef.names <- strsplit(model.terms.coef.names, ",")

  # 2) get the subset of the variables that are relevant according to the formula
  # Define the variables used in the model (and hence the unique classes of partners)
  # This is typically a subset of the available variables
  # model_vars <- c("Int", unlist(unique(lapply(model.terms.coef.names, '[[', 1))))
  model_vars <- c("Int", unique(unlist(model.terms.coef.names)))
  model_vars = model_vars[model_vars %in% colnames(Xdata)]

  Xu <- unique(Xdata[,model_vars])
  Xu <- Xu[do.call(order, as.data.frame(Xu)),]
  Zu <- unique(Zdata[,model_vars])
  Zu <- Zu[do.call(order, as.data.frame(Zu)),]
  
  Xtype <- rep(NA,nrow(Xdata))
  for(i in 1:nrow(Xu)){
    Xtype[apply(Xdata[,model_vars], 1, function(x) identical(x, Xu[i,]))] <- i
  }
  # Ztype: group membership for men (one for each man in the pop)
  Ztype <- rep(NA,nrow(Zdata))
  for(i in 1:nrow(Zu)){
    Ztype[apply(Zdata[,model_vars], 1, function(x) identical(x, Zu[i,]))] <- i
  }
  
  # order the data by pair
  # Xtype_paired = Xtype[unlist(apply(mu, 2, function(x) which(x>0)))] # for regular matrix
  Xtype_paired = Xtype[mu@i+1] # for sparse matrix
  Ztype_paired = Ztype[as.logical(colSums(mu))]
  
  # Xtype_single = table(Xtype[!rowSums(mu)])
  # Ztype_single = table(Ztype[!colSums(mu)])
  Xtype_single = table(factor(Xtype[!rowSums(mu)], 1:nrow(Xu))) # account for missing types
  Ztype_single = table(factor(Ztype[!colSums(mu)], 1:nrow(Zu))) # account for missing types
  
  pmfW=table(Xtype)/nrow(mu)
  pmfM=table(Ztype)/ncol(mu)
  
  num_Xu = nrow(Xu)
  num_Zu = nrow(Zu)
  pmfj = matrix(0,nrow=1+num_Xu, ncol=1+num_Zu) # women (X) indexed by row, men (Z) indexed by column
  pmfj[1:num_Xu,1:num_Zu] = unclass(table(Xtype_paired,Ztype_paired)) *2
  if (length(Xtype_single) > 0) {
    pmfj[1:num_Xu,1+num_Zu] = Xtype_single
  } 
  if (length(Ztype_single) > 0) {
    pmfj[1+num_Xu,1:num_Zu] = Ztype_single
  }
  
  
  pmfj <- pmfj / (nrow(Xdata) + nrow(Zdata)) # original code only divide by nrow(Xdata)?
  
  modelmat <- rpm.model.matrix(model.terms.names, model.terms.coef.names, Xu, Zu)
  
  # assume same # of explanatory variables for both men's and women's side
  NumGammaW <- num_Xu
  NumGammaM <- num_Zu
  NumGamma <- NumGammaW + NumGammaM
  NumBeta = length(theta) - NumGamma
  beta <- theta[1:NumBeta]
  GammaW <- theta[(NumBeta+1):(NumBeta+NumGammaW)]
  GammaM <- theta[(NumBeta+NumGammaW+1):length(theta)]
  
  # get the proportion of men and women
  gw = log(nrow(Xdata)/n)
  gm = log(nrow(Zdata)/n)
  
  CP_result=logitProb_latent_opp_set(beta, GammaW, GammaM, modelmat$X, modelmat$Z, gw, gm, n, symmetric)
  
  CP_result$CPW
  CP_result$CPM
  
  
  # check joint probabilities
  # couples
  pmfj_est = matrix(0,nrow=1+num_Xu, ncol=1+num_Zu)
  for (ii in 1:num_Xu) {
    for (jj in 1:num_Zu) {
      pmfj_est[ii, jj] = CP_result$CPW[ii,jj]*pmfM[jj]*exp(gm) * CP_result$CPM[jj,ii]*pmfW[ii]*exp(gw) * n
    }
  }
  # single women
  for (ii in 1:num_Xu) {
    pmfj_est[ii, 1+num_Zu] = CP_result$CPW[ii,1+num_Zu]*pmfW[ii]*exp(gw)
  }
  # single men
  for (jj in 1:num_Zu) {
    pmfj_est[1+num_Xu, jj] = CP_result$CPM[jj,1+num_Xu]*pmfM[jj]*exp(gm)
  }
  
  # debug
  # # sum of type 1 females ( ~= exp(gw) * pmfW[1] )
  # print(paste0("sum of type 1 females (should be ", exp(gw)*pmfW[1], " :"))
  # print(sum(pmfj_est[1,]))
  # 
  # # sum of type 2 females ( ~= exp(gw)  * pmfW[2] )
  # print(paste0("sum of type 2 females (should be ", exp(gw)*pmfW[2], " :"))
  # print(sum(pmfj_est[2,]))
  # 
  # # sum of all females ( ~= exp(gw) )
  # print(paste0("sum of all females (should be ", exp(gw), " :"))
  # print(sum(pmfj_est[1:2,]))
  # 
  # # sum of type 1 males ( ~= exp(gm) * pmfM[1] )
  # print(paste0("sum of type 1 males (should be ", exp(gm)*pmfM[1], " :"))
  # print(sum(pmfj_est[,1]))
  # 
  # # sum of type 2 males ( ~= exp(gm) * pmfM[2] )
  # print(paste0("sum of type 2 males (should be ", exp(gm)*pmfM[2], " :"))
  # print(sum(pmfj_est[,2]))
  # 
  # # sum of all males ( ~= exp(gm) )
  # print(paste0("sum of all males (should be ", exp(gm), " :"))
  # print(sum(pmfj_est[,1:2]))
  
  # each couple consists of 2 people
  for (ii in 1:num_Xu) {
    for (jj in 1:num_Zu) {
      pmfj_est[ii, jj] = pmfj_est[ii, jj] * 2 # CP_result$CPW[ii,jj]*pmfM[jj]*exp(gm) * CP_result$CPM[jj,ii]*pmfW[ii]*exp(gw) * n * 2
    }
  }
  
  # debug
  # # sum of population ( ~= 1 )
  # print("sum of population (should be 1")
  # print(sum(pmfj_est))
  
  # # compare estimated joint probabilities with truth
  # print("estimated joint probabilities")
  # print(pmfj_est)
  # print("observed joint probabilities")
  # print(pmfj)
  
  # # print likelihood value
  # loglik = loglikelihood_CP(beta, GammaW, GammaM, modelmat$X, modelmat$Z, pmfW, pmfM, pmfj, gw, gm, n, symmetric, sampling="INDIV")
  # logliknull = loglikelihood_null_CP(GammaW, GammaM, modelmat$X, modelmat$Z, pmfW, pmfM, pmfj, gw, gm, n, symmetric, sampling="INDIV")
  # logliknull_gamma = loglikelihood_null_CP_gamma(GammaW, GammaM, modelmat$X, modelmat$Z, pmfW, pmfM, pmfj, gw, gm, n, symmetric, sampling="INDIV")
  # print(paste0("loglik = ", loglik))
  # print(paste0("logliknull = ", logliknull))
  # # print(paste0("logliknull_gamma = ", logliknull_gamma))
  
  # # chi-squared statistics
  # if (symmetric) {
  #   df = NumBeta
  # } else {
  #   df = NumBeta # same number of explanatory variables for both sides
  # }
  # # p.val1 <- 1 - pchisq(-2*(loglik-logliknull), df=df) 
  # # print(p.val1)
  # p.val2 <- pchisq(-2*(loglik-logliknull), df=df, lower.tail=FALSE) 
  # print(paste0("chi-squared p-val = ", p.val2))
  
  # # debug
  # diffVecnull_gamma = equality_constraint_null_CP_gamma(GammaW, GammaM, modelmat$X, modelmat$Z, pmfW, pmfM, gw, gm, n, symmetric, sampling="INDIV")
  # diffVecnull = equality_constraint_null_CP(GammaW, GammaM, modelmat$X, modelmat$Z, pmfW, pmfM, gw, gm, n, symmetric, sampling="INDIV")
  # print("diffVecnull_gamma ")
  # print(diffVecnull_gamma)
  # print("diffVecnull")
  # print(diffVecnull)
  
  return(list(pmfj_est=pmfj_est, pmfj_obs=pmfj))
  
}

find_best_theta = function(beta_med, beta_sim, numBeta, numGamma) {
  best_idx = 0
  best_diff = .Machine$double.xmax
  for (ii in 1:nrow(beta_sim)) {
    diffval = sum((beta_sim[ii,1:numBeta] - beta_med)^2)
    if (diffval < best_diff) {
      best_diff = diffval
      best_idx = ii
    }
  }
  
  print(paste0("best theta [", best_idx, "] = "))
  print(beta_sim[best_idx,1:(numBeta+numGamma)])
  
  return(beta_sim[best_idx,1:(numBeta+numGamma)])
}