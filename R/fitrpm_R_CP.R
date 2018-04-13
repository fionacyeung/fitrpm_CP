#' Estimate the parameters of a Revealed Preference Matchings Model
#' 
#' \code{\link{fitrpm}} estimates parameters for men and women of certain
#' characteristics (or shared characteristics) of people of the opposite sex.
#' It does this using an approximate likelihood based on ideas from Menzel (2015).
#' 
#' @aliases rpm.object
#' @param ff formula; an \code{\link{formula}} object, of the form \code{
#' ~ <model terms>}. For the details on the possible \code{<model terms>}, see
#' \code{\link{rpm-terms}}.
#' @param mu The observed matching matrix, where 1 represent a pairing, 0 otherwise. 
#' Each row is a woman, each column is a man. The order of the rows (columns)
#' needs to be the same as in \code{Xdata} (\code{Zdata}). 
#' @param Xdata Feature matrix for women. Each row is a woman, each column is a feature. 
#' The number of column is assumed to be the same as in \code{Zdata}. 
#' @param Zdata Feature matrix for men. Each row is a man, each column is a feature. 
#' The number of column is assumed to be the same as in \code{Xdata}.
#' @param theta_0 vector; initial parameter values to initiate the search for the MLE. 
#' The length of the vector corresponds to the number of explanatory variables (number 
#' of columns in \code{Xdata}). Default vaule is a vector where the utility coefficients
#' are set to 0 and the expected utility of the opportunity set (for each distinct type of 
#' women and men) set to 1.
#' @param control A list of control parameters for algorithm tuning. Constructed using
#' \code{\link{control.fitmle}}. The \code{symmetric} parameter, if set to \code{TRUE}, indicates 
#' that the same utility coefficients are to be used for both the women's and men's sides; 
#' otherwise, they will be estimated separately for each side. The \code{sampling_protocol} parameter 
#' take the following values: \code{"INDV"} (default) (individuals are sampled, data contains both
#' singles and couples), \code{"COUPLE"} (only couples are included in the data), 
#' \code{"HOUSEHOLD"} (households are sampled, each household can be a single or a couple)
#' @return \code{\link{fitrpm}} returns an object of class \code{\link{rpm.object}}
#' that is a list consisting of the following elements: 
#' \item{coef}{The maximum likelihood estimate of \eqn{\theta}, the vector of
#' coefficients for the model parameters. This includes the model \eqn{\beta} and the model \eqn{\Gamma}.}
#' \item{loglik}{The value of the maximized log-likelihood.}
#' \item{exitflag}{integer value with the status of the optimization (4 is success).}
#' \item{call}{the call that was made to nloptr.}
#' \item{x0}{vector with starting values for the optimization.}
#' \item{message}{more informative message with the status of the optimization.}
#' \item{iterations}{number of iterations that were executed.}
#' \item{objective}{value if the objective function in the solution.}
#' \item{solution}{optimal value of the controls.}
#' \item{version}{version of NLopt that was used.}
#' \item{covar}{Approximate covriance matrix of the estimates.}
#' \item{eq}{Values from the equality constraints. Larger values indicate non-convergence.}
#' @seealso control.fitrpm, summary.rpm, print.rpm
#' @references Menzel, Konrad (2015).
#' \emph{Large Matching Markets as Two-Sided Demand Systems}
#' Econometrica, Vol. 83, No. 3 (May, 2015), 897-941.
#' @keywords models
#' @examples

library(abind)
library(nloptr)
library(numDeriv)
library(MASS)

fitrpm_R_CP <- function(formula, mu, Xdata, Zdata, theta_0=NULL, control){

    symmetric = control[["symmetric"]]
    sampling = control[["sampling_protocol"]]
    
    # number of pairs in the data
    K=nrow(Xdata)
    
    # get the proportion of men and women
    n = nrow(Xdata) + nrow(Zdata)
    gw = log(nrow(Xdata)/n)
    gm = log(nrow(Zdata)/n)
    
    # 1) parse the formula
    # intercept is always added as the first column
    Xdata <- cbind(1, Xdata)
    Zdata <- cbind(1, Zdata)
    colnames(Xdata)[1] <- "Int"
    colnames(Zdata)[1] <- "Int"
    
    model.terms <- rownames(attr(terms.formula(formula), "factors"))
    temp <- strsplit(model.terms, "[(]")
    model.terms.names <- unlist(lapply(temp, `[[`, 1))
    temp2 <- gsub(")", "", gsub("\"", "", unlist(temp)))
    temp2 <- gsub("[ \t\n\r\f\v]", "", temp2)
    model.terms.coef.names <- temp2[-which(temp2 %in% model.terms.names)]
    model.terms.coef.names <- strsplit(model.terms.coef.names, ",")
    
    # 2) get the subset of the variables that are relevant according to the formula
    # Define the variables used in the model (and hence the unique classes of partners)
    # This is typically a subset of the available variables
    model_vars <- c("Int", unlist(unique(lapply(model.terms.coef.names, '[[', 1))))
    
    
    # 3) Compute the marginal distributions of women's and men's types
    Xu <- unique(Xdata[,model_vars])
    Xu <- Xu[do.call(order, as.data.frame(Xu)),]
    Zu <- unique(Zdata[,model_vars])
    Zu <- Zu[do.call(order, as.data.frame(Zu)),]
    
    
    # 4) Create joint PMF
    # Xtype: group membership for women (one for each woman in the pop)
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
    Xtype_paired = Xtype[unlist(apply(mu, 2, function(x) which(x>0)))]
    Ztype_paired = Ztype[as.logical(colSums(mu))]
    
    Xtype_single = table(Xtype[!rowSums(mu)])
    Ztype_single = table(Ztype[!colSums(mu)])
    
    pmfW=table(Xtype)/nrow(mu)
    pmfM=table(Ztype)/ncol(mu)
    
    num_Xu = nrow(Xu)
    num_Zu = nrow(Zu)
    
    if (sampling == "COUPLE") { 
      
      pmfj = matrix(0,nrow=num_Xu, ncol=num_Zu) # women (X) indexed by row, men (Z) indexed by column
      pmfj = unclass(table(Xtype_paired,Ztype_paired))
      pmfj = pmfj / nrow(Xdata)
      
      
    } else if (sampling == "HOUSEHOLD") {
      
      pmfj = matrix(0,nrow=1+num_Xu, ncol=1+num_Zu) # women (X) indexed by row, men (Z) indexed by column
      pmfj[1:num_Xu,1:num_Zu] = unclass(table(Xtype_paired,Ztype_paired))
      
      if (length(Xtype_single) > 0) {
        pmfj[1:num_Xu,1+num_Zu] = Xtype_single
      } 
      if (length(Ztype_single) > 0) {
        pmfj[1+num_Xu,1:num_Zu] = Ztype_single
      }
      
      pmfj = pmfj / nrow(Xdata)
      
      
    } else { # assume "INDIV"
      
      pmfj = matrix(0,nrow=1+num_Xu, ncol=1+num_Zu) # women (X) indexed by row, men (Z) indexed by column
      pmfj[1:num_Xu,1:num_Zu] = unclass(table(Xtype_paired,Ztype_paired)) *2
      
      if (length(Xtype_single) > 0) {
        pmfj[1:num_Xu,1+num_Zu] = Xtype_single
      } 
      if (length(Ztype_single) > 0) {
        pmfj[1+num_Xu,1:num_Zu] = Ztype_single
      }
      
      pmfj = pmfj / (nrow(Xdata) + nrow(Zdata)) 
      
    }
   

    
    # 5) create model matrix
    modelmat <- rpm.model.matrix(model.terms.names, model.terms.coef.names, Xu, Zu)
    
    X <- modelmat$X
    Z <- modelmat$Z
    
  

    NumBetaW <- dim(X)[3]
    NumBetaM <- if(symmetric)0 else dim(Z)[3]
    NumBeta <- NumBetaW + NumBetaM
    NumGammaW <- num_Xu
    NumGammaM <- num_Zu
    NumGamma <- NumGammaW + NumGammaM
    
    
    # 6) Set theta_0 to a good starting value to save time
    if(is.null(theta_0)){
        theta_0 <- rep(c(0,1),c(NumBeta,NumGamma))
    }
    if(symmetric) {
      nstr = gsub("b1","", modelmat$Xnames)
    }else {
      nstr = c(modelmat$Xnames, modelmat$Znames)
    }
    names(theta_0) <- c(nstr, paste("ExpUtil.b1.",(1:NumGammaW),sep=""), paste("ExpUtil.b2.",(1:NumGammaM),sep="") )
    
    # core algorithm begins
    
    betaW <- theta_0[1:NumBetaW]
    betaM <- if(symmetric)NULL else theta_0[NumBetaW + (1:NumBetaM)]
    
    
    
    loglikfun <- function(theta, Xd, Zd, NumGammaW, pmfW, pmfM, pmfj, gw, gm, n, symmetric, sampling){
      if (symmetric) NumBeta <- dim(Xd)[3] else NumBeta <- dim(Xd)[3]+dim(Zd)[3]
      beta <- theta[1:NumBeta]
      GammaW <- theta[(NumBeta+1):(NumBeta+NumGammaW)]
      GammaM <- theta[(NumBeta+NumGammaW+1):length(theta)]
      -loglikelihood_CP(beta, GammaW, GammaM, Xd, Zd, pmfW, pmfM, pmfj, gw, gm, n, symmetric, sampling)
    } 
    eqfun <- function(theta, Xd, Zd, NumGammaW, pmfW, pmfM, pmfj, gw, gm, n, symmetric, sampling){
      if (symmetric) NumBeta <- dim(Xd)[3] else NumBeta <- dim(Xd)[3]+dim(Zd)[3]
      beta <- theta[1:NumBeta]
      GammaW <- theta[(NumBeta+1):(NumBeta+NumGammaW)]
      GammaM <- theta[(NumBeta+NumGammaW+1):length(theta)]
      equality_constraint_CP(beta, GammaW, GammaM, Xd, Zd, pmfW, pmfM, gw, gm, n, symmetric, sampling)
    }
    gloglikfun <- function(theta, Xd, Zd, NumGammaW, pmfW, pmfM, pmfj, gw, gm, n, symmetric, sampling){
      nl.grad(theta, loglikfun,Xd=Xd,Zd=Zd,NumGammaW=NumGammaW,
              pmfW=pmfW,pmfM=pmfM,pmfj=pmfj,gw=gw, gm=gm, n=n, symmetric=symmetric, sampling=sampling)
    }
    jeqfun <- function(theta, Xd, Zd, NumGammaW, pmfW, pmfM, pmfj, gw, gm, n, symmetric, sampling){
      nl.jacobian(theta, eqfun,Xd=Xd,Zd=Zd,NumGammaW=NumGammaW,
                  pmfW=pmfW,pmfM=pmfM, pmfj=pmfj,gw=gw,gm=gm,n=n,symmetric=symmetric, sampling=sampling)
    }
    
    # NULL likelihood
    loglikNULLfun <- function(theta, Xd, Zd, NumGammaW, pmfW, pmfM, pmfj, gw, gm, n, symmetric, sampling){
      NumBeta <- 0
      GammaW <- theta[(NumBeta+1):(NumBeta+NumGammaW)]
      GammaM <- theta[(NumBeta+NumGammaW+1):length(theta)]
      -loglikelihood_null_CP(GammaW, GammaM, Xd, Zd, pmfW, pmfM, pmfj, gw, gm, n, symmetric, sampling)
    }
    eqNULLfun <- function(theta, Xd, Zd, NumGammaW, pmfW, pmfM, pmfj, gw, gm, n, symmetric, sampling){
      NumBeta <- 0
      GammaW <- theta[(NumBeta+1):(NumBeta+NumGammaW)]
      GammaM <- theta[(NumBeta+NumGammaW+1):length(theta)]
      equality_constraint_null_CP(GammaW, GammaM, Xd, Zd, pmfW, pmfM, gw, gm, n, symmetric, sampling)
    }
    gloglikNULLfun <- function(theta, Xd, Zd, NumGammaW, pmfW, pmfM, pmfj, gw, gm, n, symmetric, sampling){
      nl.grad(theta, loglikNULLfun,Xd=Xd,Zd=Zd,NumGammaW=NumGammaW,
              pmfW=pmfW,pmfM=pmfM,pmfj=pmfj,gw=gw, gm=gm, n=n, symmetric=symmetric, sampling=sampling)
    }
    jeqNULLfun <- function(theta, Xd, Zd, NumGammaW, pmfW, pmfM, pmfj, gw, gm, n, symmetric, sampling){
      nl.jacobian(theta, eqNULLfun, Xd=Xd,Zd=Zd,NumGammaW=NumGammaW,
                  pmfW=pmfW,pmfM=pmfM, pmfj=pmfj,gw=gw,gm=gm,n=n,symmetric=symmetric, sampling=sampling)
    }
    
    
    # for C++ implementation)
  
    # eqfun <- cmpfun(eqfun)
    # loglikfun <- cmpfun(loglikfun)
    # gloglikfun <- cmpfun(gloglikfun)
    # jeqfun <- cmpfun(jeqfun)
    # 
    # eqNULLfun <- cmpfun(eqNULLfun)
    # loglikNULLfun <- cmpfun(loglikNULLfun)
    # gloglikNULLfun <- cmpfun(gloglikNULLfun)
    # jeqNULLfun <- cmpfun(jeqNULLfun)
    
    
    if(control[["algorithm"]]!="solnp"){
        out <- nloptr(x0=theta_0, eval_f=loglikfun, eval_grad_f=gloglikfun,
                      eval_g_eq=eqfun, eval_jac_g_eq=jeqfun,
                      lb=rep(c(-Inf,0),c(NumBeta,NumGamma)), # ub=rep(c(Inf, 1.0e-6),c(NumBeta,NumGamma)),
                      Xd=X,Zd=Z,NumGammaW=NumGammaW,
                      pmfW=pmfW, pmfM=pmfM, pmfj=pmfj, gw=gw, gm=gm, n=n, symmetric=symmetric, sampling=sampling,
                      opts=control)
        names(out$solution) <- names(theta_0)
        th_hat <- out$solution
        if(!is.null(names(theta_0))){names(th_hat) <- names(theta_0)}
        out$loglik <- -K*out$objective
        out$exitflag <- out$status
        cat("eq values:\n")
        print((eqfun(out$solution, X, Z, NumGammaW,pmfW, pmfM, pmfj, gw, gm, n, symmetric, sampling)))
        print(round(th_hat,2))
        
    }else{
        out <- solnp(pars=theta_0, fun=loglikfun, eqfun=eqfun,
                     LB=rep(c(-Inf,0),c(NumBeta,NumGamma)), UB=rep(c(Inf,3000),c(NumBeta,NumGamma)),
                     mu=mu,Xd=X,Zd=Z,NumGammaW=NumGammaW,
                     pmfj=pmfj,pmfW=pmfW, pmfM=pmfM,
                     control=list(tol=control[["xtol_rel"]],trace=control[["print_level"]],
                                  outer.iter=control[["maxeval"]]))
        names(out$pars) <- names(theta_0)
        th_hat <- out$pars
        out$solution <- out$pars
        #        names(th_hat) <- c(paste("betaW",(1:NumBeta)-1,sep=""), 
        #                           paste("GammaW",(1:NumGammaW)-1,sep=""), paste("GammaM",(1:NumGammaM)-1,sep="") )
        if(!is.null(names(theta_0))){names(th_hat) <- names(theta_0)}
        out$loglik <- -K*out$values[length(out$values)]
        out$exitflag <- out$converged
        cat("eq values:\n")
        print((eqfun(out$solution, X, Z, NumGammaW,pmfW, pmfM, pmfj, gw, gm, n, symmetric, sampling)))
        print(round(th_hat,2))
        
    }

    out$eq = eqfun(out$solution, X, Z, NumGammaW,pmfW, pmfM, pmfj, gw, gm, n, symmetric, sampling)
    out$coef <- th_hat[1:NumBeta]
    out$rPMFW <- 1/ (1+th_hat[(NumBeta+1):(NumBeta+NumGammaW)])
    out$rPMFM <- 1/ (1+th_hat[(NumBeta+NumGammaW+1):(NumBeta+NumGammaW+NumGammaM)])

    out$loglik <- K*loglikelihood_CP(th_hat[1:NumBeta],GammaW=th_hat[(NumBeta+1):(NumBeta+NumGammaW)], GammaM=th_hat[(NumBeta+NumGammaW+1):(NumBeta+NumGamma)],
                                     Xd=X, Zd=Z, pmfW=pmfW, pmfM=pmfM, pmfj=pmfj, gw=gw, gm=gm, n=n, symmetric=symmetric, sampling=sampling)

    
    # null likelihood
    out.null <- nloptr(x0=out$solution[-c(1:NumBeta)], eval_f=loglikNULLfun, eval_grad_f=gloglikNULLfun,
                       eval_g_eq=eqNULLfun, eval_jac_g_eq=jeqNULLfun,
                       lb=rep(0,NumGamma), # ub=rep(Inf,NumGamma),
                       Xd=X,Zd=Z,NumGammaW=NumGammaW,
                       pmfW=pmfW, pmfM=pmfM, pmfj=pmfj, gw=gw, gm=gm, n=n, symmetric=symmetric, sampling=sampling,
                       opts=control)
    out$loglik.null <- -K*out.null$objective
    #       out$loglik <- out$loglik - out$loglik.null + K*loglik.ref
    #       out$loglik.null <- K*loglik.ref
    out$pmfj=pmfj; out$pmfW=pmfW; out$pmfM=pmfM
    out$Xd=X; out$Zd=Z
    
    if(control[["hessian"]]){
      
      H <- K*hessian(loglikelihood_CP,th_hat[1:NumBeta],
                     GammaW=th_hat[(NumBeta+1):(NumBeta+NumGammaW)],GammaM=th_hat[(NumBeta+NumGammaW+1):(NumBeta+NumGamma)],
                     Xd=X,Zd=Z,pmfW=pmfW, pmfM=pmfM,pmfj=pmfj, gw=gw, gm=gm, n=n, symmetric=symmetric, sampling=sampling)
      
      Geqfun <- K*nl.jacobian(th_hat[1:NumBeta], equality_constraint_CP,  
                              GammaW=th_hat[(NumBeta+1):(NumBeta+NumGammaW)],GammaM=th_hat[(NumBeta+NumGammaW+1):(NumBeta+NumGamma)],
                              Xd=X,Zd=Z, pmfW=pmfW, pmfM=pmfM, gw=gw, gm=gm, n=n, symmetric=symmetric, sampling=sampling)
      
      dimnames(H) <- list(names(th_hat[1:NumBeta]),names(th_hat[1:NumBeta]))
      dimnames(Geqfun) <- list(names(th_hat)[(NumBeta+1):(NumBeta+NumGamma)],names(th_hat[1:NumBeta]))
      Hi <- try(ginv(-H))
      if(inherits(Hi,"try-error")){
        Hi <- -diag(1/diag(H))
      }
      V <- try(Hi - Hi %*% t(Geqfun) %*% ginv(Geqfun %*% Hi %*% t(Geqfun))%*% Geqfun %*% Hi)
      if(inherits(V,"try-error")){
        V <- Hi
      }
      if(all(is.na(diag(V)) | abs(diag(V))<1e-8)){
        out$covar <- Hi
      }else{
        out$covar <- V
      }
    }else{
      out$covar <- diag(rep(NA,NumBeta))
    }
    
    out$control <- control
    out$df <- K
    
    out$NumBeta <- NumBeta
    out$NumBetaW <- NumBetaW
    out$NumBetaM <- NumBetaM
    out$NumGammaW <- NumGammaW
    out$NumGammaM <- NumGammaM
    out$NumGamma <- NumGamma
    
    class(out) <- "rpm"
    
    return(out)
}
