# This file is the test driver. It calls the entry function (fitrpm_R_CP) for the two-sided market estimation algo.

# source("one_to_one_matching_CP.R")

library(tictoc)
library(ggplot2)

source("Gale_Shapley.R")
source("fitrpm_R_CP.R")
source("loglikelihood_CP.R")
source("equality_constraint_CP.R")
source("rpm.model.matrix.R")
source("choice_probability.R")
source("check_CP.R")

set.seed(1234)

reme = function(X,m,n){
  ##R equivalent of repmat (matlab)
  mx = dim(X)[1]
  nx = dim(X)[2]
  matrix(t(matrix(X,mx,nx*n)),mx*m,nx*n,byrow=T)
}

################## create simulated data ##########################

# specify the formula for utilities
ff = ~ b1cov("f1") + b2cov("f1") + b1absdiff("f1",1) + b2absdiff("f1",1)
# parse the formula
model.terms <- rownames(attr(terms.formula(ff), "factors"))
temp <- strsplit(model.terms, "[(]")
model.terms.names <- unlist(lapply(temp, `[[`, 1))
temp2 <- gsub(")", "", gsub("\"", "", unlist(temp)))
temp2 <- gsub("[ \t\n\r\f\v]", "", temp2)
model.terms.coef.names <- temp2[-which(temp2 %in% model.terms.names)]
model.terms.coef.names <- strsplit(model.terms.coef.names, ",")

# symmetric beta for b1 and b2
symmetric = FALSE
# true parameter for utility generation
# beta_wt <- c(1,0,1)
# beta_mt <- c(1,0,1)
beta_wt = c(0.5,0.5,1)
beta_mt = c(0.5,0.5,1)
beta_w <- matrix(beta_wt,ncol=1)
beta_m <- matrix(beta_mt,ncol=1)

# sampling protocol
sample = "INDIV" # sampling individuals
# sample = "COUPLE" # sampling couples only
# sample = "HOUSEHOLD" # sampling households, which can be a single individual or a couple

numGamma = 4
if (symmetric) {
  numBeta = length(beta_w)
} else {
  numBeta = length(beta_w) + length(beta_m)
}


################## algorithm parameters ###########################

control = list("algorithm"="NLOPT_LD_SLSQP", "symmetric"=symmetric, "sampling_protocol"=sample,
               "xtol_rel"=1.0e-8, "print_level"=0,"maxeval"=1000, 
               "ftol_rel"=1.0e-8,"check_derivatives"=FALSE,"ftol_abs"=1.0e-6,"hessian"=TRUE) 

# starting parameter for estimation
theta_0 = c(rep(0,numBeta), rep(1,numGamma))


############# simulation parameters and result ####################

# N <- c(20,50,100,200,500,1000,2000)
N = c(200)
B = 200

# N=c(20, 50)

remove_noncov = TRUE
# portions of men and women
# for nw=n*0.8, nm=n*0.2
# gw = -0.2231436
# gm = -1.609438

# for nw=n*0.9, nm=n*0.1
# gw = -0.1053605
# gm = -2.302585

# for nw=n, nm=n
gw = -0.6931472
gm = -0.6931472

################### store estimation results ######################

beta_mean <- matrix(0,ncol=numBeta,nrow=length(N))
beta_med <- matrix(0,ncol=numBeta,nrow=length(N))
beta_iqr <- matrix(0,ncol=numBeta,nrow=length(N))
beta_sd <- matrix(0,ncol=numBeta,nrow=length(N))
nonconvergence = matrix(0,ncol=1,nrow=length(N))
beta_all = matrix(0,ncol=numBeta, nrow=B*length(N))
Gamma_med <- matrix(0,ncol=numGamma,nrow=length(N))

###################################################################
###################################################################


tic.clearlog()
tic("n in N")
i = 1
sidx=1
eidx=sidx+B-1

for (n in N){

  beta_sim <- matrix(0,ncol=numBeta+2*numGamma+1,nrow=B)
    
  tic("b in B")
  for (b in 1:B){
    
    print(paste0("N=",n, " b=",b))

    nw = round(n*exp(gw))
    nm = round(n*exp(gm))
    
    # create X, Z data
    X = matrix(rep(0:1,times=nw/2),ncol=1)
    colnames(X) = "f1"
    Z = matrix(rep(0:1,times=nm/2),ncol=1)
    colnames(Z) = "f1"

    # create utility matrices
    # U_star = matrix(0, nrow=nw, ncol = nm)
    U_star = matrix(0, nrow=dim(X)[1], ncol = dim(Z)[1])
    # V_star = matrix(0, nrow=nm, ncol = nw)
    V_star = matrix(0, nrow=dim(Z)[1], ncol = dim(X)[1])
    modelmat <- rpm.model.matrix(model.terms.names, model.terms.coef.names, X, Z)
    for (ii in 1:dim(modelmat$X)[3]) {
      U_star = U_star + modelmat$X[,,ii] * beta_w[ii]
    }
    for (ii in 1:dim(modelmat$Z)[3]) {
      V_star = V_star + modelmat$Z[,,ii] * beta_m[ii]
    }
    # Menzel's code did t(V_star)
    V_star = t(V_star)
    
    # eta  <- -log(-log(matrix(runif(nw * nm), nw)))
    # zeta <- -log(-log(matrix(runif(nw * nm), nw)))
    eta  <- -log(-log(matrix(runif(dim(X)[1] * dim(Z)[1]), dim(X)[1])))
    zeta <- -log(-log(matrix(runif(dim(X)[1] * dim(Z)[1]), dim(X)[1])))
    
    
    # adjust for outside option (remain single)
    J <- round(sqrt(n))
    # eta0 <- reme(matrix(apply(-log(-log(matrix(runif(J * nw), J))), 2, max),ncol=1),1,nm)
    # zeta0 <- reme(matrix(apply(-log(-log(matrix(runif(J * nm), J))),2,max),nrow=1), nw, 1)
    eta0 <- reme(matrix(apply(-log(-log(matrix(runif(J * dim(X)[1]), J))), 2, max),ncol=1),1,dim(Z)[1])
    zeta0 <- reme(matrix(apply(-log(-log(matrix(runif(J * dim(Z)[1]), J))),2,max),nrow=1), dim(X)[1], 1)
    eta <- eta - eta0
    zeta <- zeta - zeta0
    
    U <- U_star + eta
    V <- V_star + zeta
    
    
    # generate the matching (W-optimal)
    # uses Menzel's GS which allows remaining single
    mu = Gale_Shapley(U,V)
    
    Xdata = X
    Zdata = Z
    # Compute MLE based on an observed matching
    out <- fitrpm_R_CP(ff, mu, Xdata, Zdata, theta_0, control=control)
    
    beta_sim[b,] <- c(out$solution, out$eq, out$loglik/out$loglik.null)
    
    # print(beta_sim[b,])
    
  }
  toc(log=TRUE, quiet=TRUE)
  
  ########################### sort out and print result ########################################
  if (remove_noncov) {
  noconv <- apply(beta_sim[,(numBeta+numGamma+1):(numBeta+numGamma+numGamma)] > 1.0e-5, 1, any)
  nonconvergence[i,1] = sum(noconv)
  } else {
    noconv = NULL
  }
  
  beta_cv <- beta_sim
  beta_cv[noconv,] <- NA
  
  beta_all[sidx:eidx,] = beta_cv[,1:numBeta]
  sidx=eidx+1
  eidx=sidx+B-1
  
  beta_med[i,] <- apply(beta_cv[,1:numBeta],2,median,na.rm=TRUE)
  beta_iqr[i,] <- 0.7413*apply(beta_cv[,1:numBeta],2,IQR,na.rm=TRUE)
  beta_mean[i,] <- apply(beta_cv[,1:numBeta],2,mean,na.rm=TRUE)
  beta_sd[i,] <- apply(beta_cv[,1:numBeta],2,sd,na.rm=TRUE)
  
  Gamma_med[i,] <- apply(beta_cv[,(1+numBeta):(numBeta+numGamma)],2,median,na.rm=TRUE)
  
  print(paste0("n=", n, ": nonconvergence ", sum(noconv), " of ", B))
  
  print("beta_sim:")
  print(beta_sim)
  print("beta_med:")
  print(beta_med)
  print("beta_iqr:")
  print(beta_iqr)
  print("beta_mean:")
  print(beta_mean)
  print("beta_sd:")
  print(beta_sd)
  print("Gamma_med:")
  print(Gamma_med)
  
  # check the reconstructed joint density for the matching generated from 
  # the estimated parameters
  
  # best_theta = find_best_theta(beta_med, beta_sim, numBeta, numGamma)
  # check_CP_latent(ff, best_theta, mu X, Z, symmetric) # used the last X and Z, mu
  pmfj = check_CP_latent(ff, beta_sim[B,1:(numBeta+numGamma)], mu, X, Z, symmetric)
  # compare estimated joint probabilities with truth
  print("estimated joint probabilities")
  print(pmfj$pmfj_est)
  print("observed joint probabilities")
  print(pmfj$pmfj_obs)
  
  i <- i+1

}
toc(log=TRUE, quiet=TRUE)

log.txt <- tic.log(format = TRUE)
log.lst <- tic.log(format = FALSE)
tic.clearlog()
timings <- unlist(lapply(log.lst, function(x) x$toc - x$tic))
writeLines(unlist(log.txt))

############################### display result ##############################
# create boxplot
beta_all_df = data.frame(Coeff = matrix(beta_all, ncol=1), 
                         beta = as.factor(rep(paste0("beta.",1:numBeta), each=B*length(N))), 
                         N = as.factor(rep(N, each=B)))

# excluding runs with coeff values more than 1 SD away from mean
temp2 = NULL
for (n in N) {
  for (ii in 1:numBeta) {
    idx = which(beta_all_df$beta == paste0("beta.",ii) & beta_all_df$N == n)
    temp_df = beta_all_df[idx,]
    
    betaMean = mean(temp_df$Coeff, na.rm = T)
    betaSd = sd(temp_df$Coeff, na.rm = T)
    # only extract the runs within 2 standard deviations
    beta_sub_df = subset(temp_df, Coeff < betaMean + 1*betaSd)
    beta_sub_df = subset(beta_sub_df, Coeff > betaMean + (-1)*betaSd)
    
    temp2 = rbind(temp2, beta_sub_df)
  }
}
beta_all_df = temp2

if (symmetric) {
  print(ggplot(data=subset(beta_all_df, !is.na(Coeff)), aes(x=beta, y=Coeff, fill=N)) + geom_boxplot(position=position_dodge(1)) + ggtitle("trimmed 1 SD"))
  
  png(paste0("beta_estimates_B" ,B, "_symmetric_trimmed_1SD.png"))
  print(ggplot(data=subset(beta_all_df, !is.na(Coeff)), aes(x=beta, y=Coeff, fill=N)) + geom_boxplot(position=position_dodge(1)) + ggtitle("trimmed 1 SD"))
  dev.off()
  
  beta_all_df_sub = subset(beta_all_df, Coeff < 5)
  beta_all_df_sub = subset(beta_all_df_sub, Coeff > -5)
  # see http://r.789695.n4.nabble.com/ggplot2-facet-grid-with-different-vertical-lines-on-each-facet-td3034783.html and
  # https://stackoverflow.com/questions/11846295/how-to-add-different-lines-for-facets
  hline.dat = data.frame(beta=paste0("beta.", 1:3), Coeff=beta_wt, agent.type=rep("female",3))
  
  print(ggplot(data=subset(beta_all_df_sub, !is.na(Coeff)), aes(x=N, y=Coeff)) + geom_boxplot() + facet_grid(agent.type ~ beta) +
          geom_hline(data = hline.dat, aes(yintercept =Coeff), color="red") + ggtitle("truncated (5)"))
  
  png(paste0("beta_estimates_B" ,B, "_symmetric_truncated_5.png"))
  print(ggplot(data=subset(beta_all_df_sub, !is.na(Coeff)), aes(x=N, y=Coeff)) + geom_boxplot() + facet_grid(agent.type ~ beta) +
          geom_hline(data = hline.dat, aes(yintercept =Coeff), color="red") + ggtitle("truncated (5)"))
  dev.off()
  
} else {
  beta_all_df$agent.type = "female"
  male_idx = (beta_all_df$beta %in% paste0("beta.",4:6))
  beta_all_df$agent.type[male_idx] = "male"
  
  print(ggplot(data=subset(beta_all_df, !is.na(Coeff)), aes(x=beta, y=Coeff, fill=N)) + geom_boxplot(position=position_dodge(1)))
  
  beta_all_df$beta[beta_all_df$beta == "beta.4"] = "beta.1"
  beta_all_df$beta[beta_all_df$beta == "beta.5"] = "beta.2"
  beta_all_df$beta[beta_all_df$beta == "beta.6"] = "beta.3"
  # see http://r.789695.n4.nabble.com/ggplot2-facet-grid-with-different-vertical-lines-on-each-facet-td3034783.html and
  # https://stackoverflow.com/questions/11846295/how-to-add-different-lines-for-facets
  hline.dat = data.frame(beta=paste0("beta.",rep(1:3,times=2)), Coeff=beta_wt, agent.type=c(rep("female",3), rep("male",3)))
  
  
  beta_all_df_sub = subset(beta_all_df, Coeff < 3)
  beta_all_df_sub = subset(beta_all_df_sub, Coeff > -3)
  print(ggplot(data=subset(beta_all_df_sub, !is.na(Coeff)), aes(x=N, y=Coeff)) + geom_boxplot() + facet_grid(agent.type ~ beta) +
    geom_hline(data = hline.dat, aes(yintercept =Coeff), color="red") + ggtitle("truncated (3)"))
  
  png(paste0("B",B,"_beta_estimates_truncated_3.png"))
  print(ggplot(data=subset(beta_all_df_sub, !is.na(Coeff)), aes(x=N, y=Coeff)) + geom_boxplot() + facet_grid(agent.type ~ beta) +
          geom_hline(data = hline.dat, aes(yintercept =Coeff), color="red") + ggtitle("truncated (3)"))
  dev.off()
  
  print(ggplot(data=subset(beta_all_df, !is.na(Coeff)), aes(x=N, y=Coeff)) + geom_boxplot() + facet_grid(agent.type ~ beta) +
    geom_hline(data = hline.dat, aes(yintercept =Coeff), color="red") + ggtitle("trimmed 1 SD"))
  
  png(paste0("B",B,"_beta_estimates_trimmed_1_SD.png"))
  print(ggplot(data=subset(beta_all_df, !is.na(Coeff)), aes(x=N, y=Coeff)) + geom_boxplot() + facet_grid(agent.type ~ beta) +
          geom_hline(data = hline.dat, aes(yintercept =Coeff), color="red") + ggtitle("trimmed 1 SD"))
  dev.off()
}




